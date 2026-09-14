# Paired-end bulk RNA-seq Nextflow pipeline

A compact teaching/research scaffold using Nextflow DSL2.

## Workflow

1. Raw-read QC: FastQC
2. Adapter and quality trimming: fastp
3. Genome index and alignment: STAR
4. Gene-level counting: featureCounts
5. Two-condition differential expression: DESeq2
6. Combined QC report: MultiQC

## Requirements

- Linux or macOS
- Java 17 or later
- Nextflow 24.10 or later
- Conda or Mamba

## Input

Edit `assets/samplesheet.csv`. Required columns:

- `sample`: unique sample name, matching the desired count-matrix column
- `condition`: exactly two groups for the included DESeq2 template
- `fastq_1`, `fastq_2`: paths to paired-end gzipped FASTQ files

Use a genome FASTA and a matching GTF annotation from the same reference build and provider.

## Run

```bash
nextflow run main.nf \
  --input assets/samplesheet.csv \
  --fasta reference/genome.fa \
  --gtf reference/annotation.gtf \
  --outdir results \
  -profile conda \
  -resume
```

Skip differential expression when you only need QC, alignment and counts:

```bash
nextflow run main.nf --input assets/samplesheet.csv --fasta reference/genome.fa --gtf reference/annotation.gtf --skip_deseq2 true -profile conda
```

For SLURM, adjust resources in `nextflow.config`, then use `-profile slurm`.

## Main outputs

- `results/qc/raw_fastqc/`
- `results/trimmed/`
- `results/alignment/`
- `results/counts/`
- `results/deseq2/`
- `results/multiqc/multiqc_report.html`
- Nextflow execution report, timeline, trace and DAG in the launch directory

## Important assumptions

- Paired-end bulk RNA-seq only.
- `featureCounts` uses `-p --countReadPairs`. Confirm strandedness and add `-s 1` or `-s 2` when appropriate.
- DESeq2 compares the second factor level with the first. For explicit reference selection, edit `bin/deseq2.R` and use `relevel()`.
- Biological replicates are required for meaningful differential-expression inference.
- Run first on a small test dataset before scaling.
