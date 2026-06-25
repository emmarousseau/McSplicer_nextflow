# McSplicer Nextflow Pipeline

Nextflow DSL2 pipeline for splice site usage estimation using the McSplicer probabilistic model.

## Quick Start

```bash
# Local execution
nextflow run main.nf --input samplesheet.tsv --outdir results

# With Docker
nextflow run main.nf --input samplesheet.tsv --outdir results -profile docker

# On SLURM cluster
nextflow run main.nf --input samplesheet.tsv --outdir results -profile slurm

# Resume after failure
nextflow run main.nf --input samplesheet.tsv --outdir results -resume

# Test with bundled example data
nextflow run main.nf -profile test
```

## Pipeline Steps

1. **exonRefine** — Refines GTF annotation into non-overlapping subexonic regions
2. **sigcount** — Counts reads mapped to signatures and junctions (`.cnt` file)
3. **McSplicer** — EM algorithm estimating splice site usage probabilities

## Samplesheet Format

Tab-separated file with the following columns:

| Column | Required | Description |
|--------|----------|-------------|
| `sample_id` | yes | Unique sample identifier |
| `bam_file` | yes | Path to BAM file (must be indexed) |
| `gtf_file` | yes | Path to GTF annotation file |
| `read_len` | yes | Read length (integer) |
| `gene_id` | no | Restrict to a specific gene (e.g., ENSG00000145293) |
| `prefix` | no | Output file prefix/suffix |
| `paired_end` | no | Use paired-end read counting (true/false, default: false) |
| `skip_preprocess` | no | Skip exonRefine+sigcount (true/false, default: false) |
| `bootstraps` | no | Number of bootstrap replicates (default: 0) |
| `anno` | no | Include transcript annotation in output (y/n, default: n) |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | — | Path to samplesheet TSV (required) |
| `--outdir` | `./results` | Output directory |
| `--read_len` | 151 | Default read length |
| `--bootstraps` | 0 | Default bootstrap replicates |
| `--anno` | n | Include transcript-to-splice-site mapping |
| `--paired_end` | false | Default paired-end mode |
| `--skip_preprocess` | false | Skip preprocessing globally |

## Profiles

| Profile | Description |
|---------|-------------|
| `standard` | Local execution (default) |
| `docker` | Run in Docker container |
| `singularity` | Run in Singularity container |
| `conda` | Use conda environment |
| `slurm` | Submit to SLURM scheduler |
| `lsf` | Submit to LSF scheduler |
| `test` | Run with bundled example data |

## Output Structure

```
results/
├── refined_gtf/             # Refined GTF files from exonRefine
├── signature_counts/        # Read count files (.cnt)
├── mcsplicer/               # Splice site usage CSV files
│   └── output/{chr}/{gene}.csv
└── pipeline_info/           # Execution reports and trace
```

## Docker

Build the container:

```bash
docker build -t mcsplicer:latest .
```

## Requirements

- Nextflow >= 23.04
- Python >= 3.6 with numpy, matplotlib
- Platform-specific binaries (exonRefine, sigcount) included in `bin/`
