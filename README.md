# PGxIntegrator

PGxIntegrator is a stage-based pipeline for preparing pharmacogenomic VCF inputs and running **Aldy4** for star-allele calling.

## Features

* Prepare platform-specific and integrated per-gene VCFs
* Split cohort-level VCFs into single-sample VCFs
* Run Aldy4 in parallel across genes and samples
* Support config-based runs with optional CLI overrides
* Organize outputs, logs, and restartable runs

## Workflow

```bash
pgx prepare --config config.yaml
pgx split --config config.yaml
pgx aldy --config config.yaml
```

You can also override config values from the command line:

```bash
pgx aldy --config config.yaml --jobs 32 --resume
```

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
  gene_manifest: /path/to/gene_manifest.tsv

prepare:
  output_dir: /path/to/prepared
  normalize: true
  harmonize_array: true

liftover:
  target_build: hg38
  enabled: true
  chain_file: /path/to/sourceToTarget.over.chain.gz

split:
  output_dir: /path/to/single_sample_vcfs
  samples_per_batch: 1000

aldy:
  output_dir: /path/to/aldy_results
  jobs: 16
  genome: hg38
  profile: wgs
  reference_fasta: /path/to/hs38DH.fa
```

### 2. Prepare per-gene VCFs

```bash
pgx prepare --config config.yaml
```

### 3. Split per-gene VCFs into single-sample VCFs

```bash
pgx split --config config.yaml
```

### 4. Run Aldy4

```bash
pgx aldy --config config.yaml
```

## Input files

### Sample manifest

Tab-delimited file listing the samples to process.

Example:

```tsv
sample_id	array_id	wes_id
S001	12345_12345	S001
S002	67890_67890	S002
```

### Gene manifest

Tab-delimited file mapping each gene and platform to an input VCF. Each input must declare its source build explicitly.

Example:

```tsv
gene	platform	build	path
CYP2B6	wes	hg38	/data/wes/CYP2B6_wes_hg38.vcf.gz
CYP2B6	array	hg19	/data/array/CYP2B6_array_hg19.vcf.gz
SLCO1B1	wes	hg19	/data/wes/SLCO1B1_wes_hg19.vcf.gz
SLCO1B1	array	hg38	/data/array/SLCO1B1_array_hg38.vcf.gz
```

If an input build differs from the pipeline target build, that input can be lifted over during `pgx prepare`.

## Commands

### `pgx prepare`

Prepares per-gene VCFs from the provided inputs.

Typical tasks may include:

* sample extraction
* sample renaming
* liftover when needed
* normalization
* harmonization
* integration of WES and array variants

Example:

```bash
pgx prepare --config config.yaml
```

With CLI overrides:

```bash
pgx prepare --config config.yaml --output-dir ./prepared
```

### `pgx split`

Splits per-gene cohort VCFs into single-sample VCFs.

Example:

```bash
pgx split --config config.yaml
```

### `pgx aldy`

Runs Aldy4 on prepared single-sample VCFs.

Example:

```bash
pgx aldy --config config.yaml
```

With CLI overrides:

```bash
pgx aldy --config config.yaml --jobs 32 --resume
```

## Output layout

Example directory structure:

```text
results/
  prepared/
    CYP2B6_merged.vcf.gz
    SLCO1B1_merged.vcf.gz
  single_sample_vcfs/
    CYP2B6/
      batch_000/
        SAMPLE1.vcf.gz
        SAMPLE1.vcf.gz.csi
    SLCO1B1/
      batch_000/
        SAMPLE1.vcf.gz
        SAMPLE1.vcf.gz.csi
  aldy_results/
    CYP2B6/
      batch_000/
        SAMPLE1_CYP2B6_aldy.out
        SAMPLE1_CYP2B6_aldy.log
```

## Notes

* Declare the source build for every input dataset.
* Use one target build consistently across all stages.
* Any input source can require liftover during `pgx prepare` if its build differs from the target build.
* Aldy settings such as `profile`, `genome`, and `reference_fasta` should match the prepared VCFs.
* Completed runs can be resumed using `--resume` when supported.

## Roadmap

Planned additions include:

* packaged release on PyPI
* example datasets
* validation checks for manifests and reference builds
* optional downstream summary utilities

## License

Add your license here.

## Citation

Citation information will be added after release.
