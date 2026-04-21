# PGxIntegrator

PGxIntegrator is a stage-based pipeline for harmonizing and integrating pharmacogenomic variant data from multiple genomic sources, including whole-exome sequencing (WES) and genotyping array data, into per-gene VCFs for downstream star-allele calling with Aldy4.

The pipeline was developed and evaluated in a UK Biobank study benchmarking WES, imputed array, and integrated call sets against whole-genome sequencing (WGS) across 33 pharmacogenes. As software, however, PGxIntegrator is intended as a reusable workflow for scalable pharmacogenomic data integration and star-allele analysis when multi-source genomic data are available.

## Features

- Harmonize and integrate PGx-relevant variants across multiple sources
- Split cohort-level per-gene VCFs into single-sample VCFs
- Run Aldy4 in parallel across genes and samples
- Support config-based runs with optional CLI overrides
- Organize outputs, logs, and restartable runs

## Workflow

```bash
pgx prepare --config config.yaml
pgx split --config config.yaml
pgx aldy --config config.yaml
```

You can also override selected config values from the command line.

Example:

```bash
pgx aldy --config config.yaml --jobs 32 --resume
```

## Dependencies

### Python

Install Python dependencies with:

```bash
pip install -r requirements.txt
```

### External tools

PGxIntegrator depends on external bioinformatics tools that must be installed separately and available on your `PATH`:

- `bcftools`
- `bgzip`
- `liftOver`
- `Aldy4`

### Versions used during development

The pipeline is not strictly limited to these versions, but the following were used during development and testing:

- Aldy4 4.6
- bcftools 1.20

## Installation

```bash
git clone <repo-url>
cd PGxIntegrator
pip install -e .
```

## Quick start

### 1. Create a config file

Example:

```yaml
project_name: example_pgx

reference:
  fasta: /path/to/hs38DH.fa
  build: hg38

inputs:
  sample_manifest: /path/to/sample_manifest.tsv
  dataset_manifest: /path/to/dataset_manifest.tsv
  gene_list: /path/to/gene_list.tsv

prepare:
  output_dir: /path/to/prepared
  normalize: true
  harmonize: true
  sample_missingness_threshold: 0.05
  keep_intermediate: false
  resume: false

liftover:
  enabled: true
  target_build: hg38
  chain_file: /path/to/sourceToTarget.over.chain.gz

split:
  output_dir: /path/to/single_sample_vcfs
  input_dir: /path/to/prepared
  samples_per_batch: 1000
  jobs: 4
  resume: false

aldy:
  output_dir: /path/to/aldy_results
  input_dir: /path/to/single_sample_vcfs
  jobs: 16
  genome: hg38
  profile: wgs
  reference_fasta: /path/to/hs38DH.fa
  resume: false
```

### 2. Prepare platform-specific and integrated per-gene VCFs

```bash
pgx prepare --config config.yaml
```

Example with optional overrides:

```bash
pgx prepare --config config.yaml --resume
```

### 3. Split per-gene cohort VCFs into single-sample VCFs

```bash
pgx split --config config.yaml
```

Example with optional overrides:

```bash
pgx split --config config.yaml --jobs 8 --resume
```

### 4. Run Aldy4

```bash
pgx aldy --config config.yaml
```

Example with optional overrides:

```bash
pgx aldy --config config.yaml --jobs 32 --resume
```

## Input files

### Sample manifest

Tab-delimited file listing the samples to process.

Example:

```tsv
sample_id	array_id	wes_id
S001	ARR01	WES01
S002	ARR02	WES02
S003	ARR03	WES03
```

### Dataset manifest

Tab-delimited file describing each input dataset.

Example:

```tsv
dataset_id	platform	build	layout	chrom	gene	path
wes_chr10	wes	hg38	cohort_by_chrom	10		../wes/chr10_wes_hg38.vcf.gz
array_chr10	array	hg19	cohort_by_chrom	10		../array/chr10_array_hg19.vcf.gz
wes_vkorc1	wes	hg38	gene_by_file		VKORC1	../wes/VKORC1_wes_hg38.vcf.gz
array_vkorc1	array	hg19	gene_by_file		VKORC1	../array/VKORC1_array_hg19.vcf.gz
```

Field descriptions:

- `dataset_id`: unique dataset name
- `platform`: input platform such as `wes` or `array`
- `build`: source genome build of the dataset
- `layout`: input layout such as `cohort_by_chrom` or `gene_by_file`
- `chrom`: chromosome for chromosome-level inputs
- `gene`: gene name for gene-level inputs
- `path`: path to the input VCF

### Gene list

Tab-delimited file listing genes to process.

Example:

```tsv
gene
CYP2C19
VKORC1
CYP2B6
SLCO1B1
```

## Commands

### `pgx prepare`

Prepare platform-specific and integrated per-gene VCFs from the provided inputs.

Typical tasks may include:

- sample extraction
- sample renaming
- liftOver when needed
- normalization
- harmonization
- integration of variants across genomic sources

Example:

```bash
pgx prepare --config config.yaml
```

### `pgx split`

Split per-gene cohort VCFs into single-sample VCFs.

Example:

```bash
pgx split --config config.yaml
```

### `pgx aldy`

Run Aldy4 on split single-sample VCFs.

Example:

```bash
pgx aldy --config config.yaml
```

Example with parameter overrides:

```bash
pgx aldy --config config.yaml --jobs 32 --resume
```

## Output layout

Example directory structure:

```text
results/
  prepared/
    CYP2C19.vcf.gz
    VKORC1.vcf.gz
  split/
    CYP2C19/
      batch_000/
        S001.vcf.gz
        S001.vcf.gz.csi
    VKORC1/
      batch_000/
        S001.vcf.gz
        S001.vcf.gz.csi
  aldy/
    CYP2C19/
      batch_000/
        S001_aldy.out
        S001_aldy.log
    VKORC1/
      batch_000/
        S001_aldy.out
        S001_aldy.log
```

## Notes

- Use one target build consistently across the pipeline.
- If an input dataset build differs from the target build, it can be lifted over during `pgx prepare`.
- Aldy settings such as `profile`, `genome`, and `reference_fasta` should match the prepared VCFs.
- Completed runs can be resumed when supported.
- Relative paths in the YAML config are resolved relative to the config file location.
- The example dataset is intended as a small release-testing dataset.

## License

MIT License