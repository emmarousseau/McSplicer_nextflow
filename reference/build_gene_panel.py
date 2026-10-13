#!/usr/bin/env python3
"""
Filter the Ensembl GRCh38 GTF down to a gene panel, for McSplicer.

Emits:
    panel_genes.gtf    - full gene/transcript/exon/... records for the panel genes,
                          fed to exonRefine.
    panel_regions.bed  - one padded region per gene (contig, start, end, gene_name),
                          fed to `samtools view -L` when subsetting BAMs.

Ensembl GTFs identify genes by `gene_id` (ENSG...), not symbol — match on the
`gene_name` attribute instead so a plain gene-symbol panel list works directly.

Usage:
    python3 build_gene_panel.py --gtf Homo_sapiens.GRCh38.116.gtf.gz --out_dir ./panel \
        [--genes_file genes.txt] [--pad 10000]
"""

import argparse
import gzip
import re
from pathlib import Path

# STMN2 + the 39-gene TDP-43 cryptic-splicing panel from
# LRL_GENI_NYGC_ALS_analysis/scripts/07_gene_exitron_correlation_table.R,
# plus TARDBP, ELAVL3, MASP2 (UNC13A/STMN2 already present) - 43 genes total.
DEFAULT_PANEL = [
    "ACTL6B", "AKT3", "ARMH1", "ATG4B", "CAMK2B", "CBARP", "CELF5",
    "CELSR3", "CHID1", "ELAVL3", "EPB41L4A", "ETS2", "FEZ1", "HDAC6", "INSR",
    "IQCE", "KALRN", "KNDC1", "LRP8", "MASP2", "MNAT1", "NSMCE4A", "NTRK3",
    "NUP188", "NYNRIN", "PFKP", "PRELID3A", "PTPRT", "PXDN", "RAB40B",
    "RAP1GAP", "SETD5", "SLC24A3", "STMN2", "TARDBP", "TMEM117", "TRAPPC12",
    "UNC13A", "UNC13B", "USP13", "WDR35", "ZNF382", "ZNF529",
]

ATTR_GENE_NAME_RE = re.compile(r'gene_name "([^"]+)"')


def open_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def load_gene_list(genes_file: Path | None) -> list[str]:
    if genes_file is None:
        return DEFAULT_PANEL
    genes = [line.strip() for line in genes_file.read_text().splitlines()]
    return [g for g in genes if g]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", type=Path, required=True, help="Reference GTF (plain or .gz)")
    parser.add_argument("--out_dir", type=Path, required=True, help="Output directory")
    parser.add_argument("--genes_file", type=Path, default=None,
                         help="Optional file, one gene symbol per line (default: built-in 43-gene panel)")
    parser.add_argument("--pad", type=int, default=10_000,
                         help="Padding (bp) added to each side of a gene's start/end for the BED (default: 10000)")
    args = parser.parse_args()

    genes = load_gene_list(args.genes_file)
    gene_set = set(genes)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    out_gtf_path = args.out_dir / "panel_genes.gtf"
    out_bed_path = args.out_dir / "panel_regions.bed"

    # contig, min_start, max_end per gene_name (1-based inclusive GTF coords)
    gene_extent: dict[str, list] = {}
    found_genes: set[str] = set()

    with open_maybe_gzip(args.gtf) as fin, open(out_gtf_path, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                continue
            m = ATTR_GENE_NAME_RE.search(line)
            if not m:
                continue
            gene_name = m.group(1)
            if gene_name not in gene_set:
                continue

            fout.write(line)
            found_genes.add(gene_name)

            fields = line.rstrip("\n").split("\t")
            contig, start, end = fields[0], int(fields[3]), int(fields[4])
            if gene_name not in gene_extent:
                gene_extent[gene_name] = [contig, start, end]
            else:
                entry = gene_extent[gene_name]
                entry[1] = min(entry[1], start)
                entry[2] = max(entry[2], end)

    with open(out_bed_path, "w") as fbed:
        for gene_name in sorted(gene_extent, key=lambda g: (gene_extent[g][0], gene_extent[g][1])):
            contig, start, end = gene_extent[gene_name]
            bed_start = max(0, start - 1 - args.pad)  # BED is 0-based half-open
            bed_end = end + args.pad
            fbed.write(f"{contig}\t{bed_start}\t{bed_end}\t{gene_name}\n")

    missing = gene_set - found_genes
    print(f"Genes requested : {len(gene_set)}")
    print(f"Genes found     : {len(found_genes)}")
    print(f"Panel GTF       : {out_gtf_path}")
    print(f"Panel BED       : {out_bed_path}")
    if missing:
        print(f"WARNING: {len(missing)} gene(s) not found in GTF: {sorted(missing)}")


if __name__ == "__main__":
    main()
