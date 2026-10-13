#!/usr/bin/env python3
"""
Stage gene-panel-subsetted BAMs for the full NYGC ALS Consortium cohort and
emit a McSplicer_nextflow samplesheet.

For each sample in coldata.csv (joined to its star_salmon BAM by SRR
accession):
  1. Presigns the S3 bam/bai (so samtools can do indexed HTTPS range reads
     without AWS creds).
  2. Runs `samtools view -b -M -X <bam_url> <bai_url> -L panel_regions.bed`
     to pull out only the gene-panel regions, producing a small local BAM.
     `-M -L` (not plain `-L`) is required over presigned HTTPS URLs — plain
     `-L` hangs indefinitely on remote URLs.
  3. Indexes the subset BAM.
  4. Reads the actual read length from the subset BAM (first alignment's
     SEQ length) rather than assuming — this cohort mixes 100bp and 125bp+
     batches.
  5. Emits a McSplicer_nextflow samplesheet.tsv row.

Usage:
    python3 stage_nygc_panel_bams.py \
        --coldata /path/to/coldata.csv \
        --s3_prefix s3://<bucket>/<prefix>/nf-core_rnaseq/<run_id> \
        --bed ../reference/panel/panel_regions.bed \
        --gtf ../reference/panel/panel_genes.gtf \
        --out_dir /path/to/scratch/nygc_panel_bams \
        --aws_profile <profile> \
        --workers 8

The 43-gene panel itself is NOT set per-sample here (the samplesheet's
gene_id column only supports a single gene) — pass it once via
--gene_list <reference/panel/panel_gene_id_list.txt> on the
`nextflow run main.nf` invocation instead.
"""

import argparse
import csv
import logging
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger(__name__)


def load_samples(coldata: Path) -> list[dict]:
    with open(coldata, newline="") as f:
        return list(csv.DictReader(f))


def list_bam_keys(s3_prefix: str, profile: str) -> dict[str, str]:
    """List star_salmon/*.markdup.sorted.bam and map SRR accession -> filename."""
    cmd = ["aws", "s3", "ls", "--profile", profile, f"{s3_prefix}/star_salmon/"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"aws s3 ls failed: {result.stderr}")

    srr_to_filename = {}
    for line in result.stdout.splitlines():
        parts = line.split()
        if not parts:
            continue
        name = parts[-1]
        if not name.endswith(".markdup.sorted.bam"):
            continue
        base = name[: -len(".markdup.sorted.bam")]
        srr = base.split(".")[-1]
        srr_to_filename[srr] = name
    return srr_to_filename


def presign(s3_path: str, profile: str, expires_in: int) -> str:
    cmd = ["aws", "s3", "presign", s3_path, "--expires-in", str(expires_in), "--profile", profile]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"aws s3 presign failed for {s3_path}: {result.stderr}")
    return result.stdout.strip()


def subset_bam(bam_url: str, bai_url: str, bed_file: Path, out_bam: Path) -> None:
    cmd = ["samtools", "view", "-b", "-M", "-X", bam_url, bai_url, "-L", str(bed_file), "-o", str(out_bam)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"samtools view failed: {result.stderr}")

    index_cmd = ["samtools", "index", str(out_bam)]
    result = subprocess.run(index_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"samtools index failed: {result.stderr}")


def read_count(bam: Path) -> int:
    result = subprocess.run(["samtools", "view", "-c", str(bam)], capture_output=True, text=True)
    if result.returncode != 0:
        return -1
    return int(result.stdout.strip())


def detect_read_len(bam: Path) -> int:
    """First aligned read's SEQ length. Falls back to 0 if the subset is empty."""
    result = subprocess.run(["samtools", "view", str(bam)], capture_output=True, text=True)
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) > 9:
            return len(fields[9])
    return 0


def process_sample(row: dict, srr_to_filename: dict, s3_prefix: str, bed_file: Path,
                    out_dir: Path, profile: str, presign_expiry: int) -> dict | None:
    srr = row["sample_id"].strip()
    filename = srr_to_filename.get(srr)
    if filename is None:
        log.error(f"[{srr}] no BAM found under star_salmon/ — skipping")
        return None

    out_bam = out_dir / f"{srr}.panel.bam"
    bam_key = f"{s3_prefix}/star_salmon/{filename}"
    bai_key = f"{bam_key}.bai"

    try:
        bam_url = presign(bam_key, profile, presign_expiry)
        bai_url = presign(bai_key, profile, presign_expiry)

        t0 = time.time()
        subset_bam(bam_url, bai_url, bed_file, out_bam)
        elapsed = time.time() - t0

        n_reads = read_count(out_bam)
        read_len = detect_read_len(out_bam)
        log.info(f"[{srr}] done in {elapsed:.1f}s, {n_reads} reads, read_len={read_len} -> {out_bam.name}")

        return {
            "sample_id": srr,
            "bam_file": str(out_bam.resolve()),
            "read_len": read_len,
            "n_reads": n_reads,
        }
    except RuntimeError as e:
        log.error(f"[{srr}] FAILED: {e}")
        return None


def write_samplesheet(rows: list[dict], gtf_file: Path, out_path: Path) -> None:
    # gene_id (singular) is left blank: the 43-gene panel is passed once via
    # params.gene_list at the `nextflow run` invocation, not per-row — the
    # per-row gene_id column only supports a single gene_id, not a list.
    #
    # prefix is set to sample_id: MCSPLICER's publishDir copies every sample's
    # output into the same flat output/{chr}/{gene_id}_{prefix}.csv path with
    # no per-sample subdirectory, so prefix must be unique per sample or
    # samples will silently overwrite each other's CSVs.
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample_id", "bam_file", "gtf_file", "read_len", "gene_id", "prefix", "paired_end", "skip_preprocess"])
        for r in rows:
            writer.writerow([
                r["sample_id"], r["bam_file"], str(gtf_file), r["read_len"], "", r["sample_id"],
                "true", "false",
            ])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coldata", type=Path, required=True, help="coldata.csv with a sample_id (SRR) column")
    parser.add_argument("--s3_prefix", type=str, required=True, help="S3 prefix of the nf-core/rnaseq run (no trailing slash)")
    parser.add_argument("--bed", type=Path, required=True, help="panel_regions.bed (from build_gene_panel.py)")
    parser.add_argument("--gtf", type=Path, required=True, help="panel_genes.gtf (from build_gene_panel.py)")
    parser.add_argument("--out_dir", type=Path, required=True, help="Output directory for subset BAMs + samplesheet.tsv")
    parser.add_argument("--aws_profile", type=str, default="local", help="AWS CLI profile for presigning")
    parser.add_argument("--presign_expiry", type=int, default=3600, help="Presigned URL expiry in seconds (default: 3600)")
    parser.add_argument("--workers", type=int, default=8, help="Parallel workers (default: 8)")
    parser.add_argument("--limit", type=int, default=None, help="Only process the first N samples (for a dry-run pilot)")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    samples = load_samples(args.coldata)
    if args.limit:
        samples = samples[: args.limit]
    log.info(f"Loaded {len(samples)} sample(s) from {args.coldata}")

    log.info("Listing star_salmon/ for BAM filenames...")
    srr_to_filename = list_bam_keys(args.s3_prefix, args.aws_profile)
    log.info(f"Found {len(srr_to_filename)} BAM(s) under star_salmon/")

    tsv_rows = []
    failed = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(process_sample, row, srr_to_filename, args.s3_prefix, args.bed,
                        args.out_dir, args.aws_profile, args.presign_expiry): row["sample_id"]
            for row in samples
        }
        for fut in as_completed(futures):
            srr = futures[fut]
            result = fut.result()
            if result is None:
                failed.append(srr)
            else:
                tsv_rows.append(result)

    samplesheet_path = args.out_dir / "samplesheet.tsv"
    write_samplesheet(tsv_rows, args.gtf.resolve(), samplesheet_path)
    log.info(f"Wrote {len(tsv_rows)} sample row(s) to {samplesheet_path}")
    if failed:
        log.warning(f"{len(failed)} sample(s) failed: {failed}")


if __name__ == "__main__":
    main()
