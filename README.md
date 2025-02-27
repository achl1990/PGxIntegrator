# PGxIntegrator

**PGxIntegrator** is a **pipeline** for **integrating Whole Exome Sequencing (WES) and imputed Array data into a WGS-like VCF** for pharmacogenomic analysis. It is a pipelining tool that integrates various bioinformatics tools for ease of use, offering **multi-sample VCF splitting, Liftover (hg19 → hg38), preprocessing, and parallel Aldy4-based star-allele genotyping**.

## Features

✅ **WES + Array Integration** → Merges WES and (imputed) array data into a single VCF.\
✅ **Multi-Sample VCF Processing** → Splits multi-sample VCFs into single samples and runs Aldy4 in parallel.\
✅ **Automated Preprocessing** → Standardizes, normalizes, and filters variants.\
✅ **Liftover Support** → Converts `hg19 → hg38` if needed.\
✅ **Flexible & Scalable** → Works with both single-sample and multi-sample VCFs. \
✅ **Fetch Gene Coordinates** → Retrieves gene coordinates from **Ensembl (GRCh37/hg19 or GRCh38/hg38)** for updating pharmacogenes.

## Installation

PGxIntegrator is a pipelining tool that automates various steps in pharmacogenomic analysis, requiring both Python packages and system tools. Python dependencies can be installed using `requirements.txt`, while system tools must be installed separately.

PGxIntegrator requires Python and several bioinformatics tools, including `bcftools`, and `Liftover`. `Aldy4` is required only if running Aldy analysis.

```bash
# Clone the repository
git clone https://github.com/achl1990/PGxIntegrator.git
cd PGxIntegrator

# Install dependencies
pip install -r requirements.txt
```

## Usage

**⚠ Note:** The `--integrate` function is currently under development and will be disabled until the next update.

### **Basic Command**

```bash
python pgxintegrator.py --wes wes.vcf --reference hs38DH.fa --gene CYP2C19  # OR
python pgxintegrator.py --array array.vcf --reference hs38DH.fa --gene CYP2C19
```

### **Full List of Arguments**

| Argument                  | Description                                             | Default                 |
| ------------------------- | ------------------------------------------------------- | ----------------------- |
| `--wes` | Path to WES VCF (or CRAM/BAM) | **Required if --array is not provided**            |
| `--array` | Path to imputed Array VCF | **Required if --wes is not provided**                |
| `--reference`             | Path to reference genome (`hs38DH.fa`)                  | **Required**            |
| `--reference-build-wes`   | Build of WES data (`hg19` or `hg38`)                    | **Required**            |
| `--reference-build-array` | Build of Array data (`hg19` or `hg38`)                  | **Required if --array** |
| `--liftover`              | Apply hg19 → hg38 LiftOver                              | `False`                 |
| `--gene`                  | Gene(s) to analyze (`CYP2C19,VKORC1, all`)              | **Required**            |
| `--integrate`             | Enable WES + Array integration (**Currently Disabled**) | `False`                 |
| `--harmonization`         | Standardize, normalize, and filter                      | `True`                  |
| `--multi-sample`          | Treat VCF as multi-sample                               | `False`                 |
| `--split-multi`           | Split multi-sample VCF                                  | `True`                  |
| `--run-aldy`              | Run Aldy analysis                                       | `True`                  |
| `--jobs`                  | Number of parallel jobs for GNU parallel                | `8`                     |
| `--aldy-output`           | Path to store Aldy results                              | `./results`             |
| `--output-dir`            | Directory for final files                               | `./results`             |
| `--log`                   | Path to log file                                        | `./log.txt`             |
| `--keep-temp`             | Retain intermediate files                               | `False`                 |
| `--jobs`                  | Number of CPU threads                                   | `8`                     |
| `--dry-run`               | Run pipeline without executing                          | `False`                 |

## Example Workflows

### **Integrate WES + Array and run Aldy**

```bash
python pgxintegrator.py --wes sample_wes.vcf --array sample_array.vcf \
  --reference hs38DH.fa --gene CYP2C19 --integrate --run-aldy
```

### **Split Multi-Sample VCF and Run Aldy in Parallel**

```bash
python pgxintegrator.py --wes multi_sample.vcf --split-multi --run-aldy --jobs 16
```

### **Fetching Gene Coordinates from Ensembl**
PGxIntegrator can retrieve genomic coordinates for pharmacogenes from **Ensembl**, helping update Aldy configurations when new genes are added.

```bash
python pgxintegrator.py --fetch-coordinates --build hg38 --output gene_coordinates.tsv
```

### **Apply Liftover from hg19 to hg38**

```bash
python pgxintegrator.py --wes sample_wes_hg19.vcf --liftover --reference hs38DH.fa
```

## Dependencies

### **Python Packages**
These can be installed using `requirements.txt`:
```bash
pip install -r requirements.txt
```

- `pandas`
- `requests`
- `matplotlib`

### **System Tools**
These must be installed separately:
```bash
# Install via apt (Ubuntu/Debian)
sudo apt install bcftools parallel

# Install via conda
conda install bioconda::bcftools bioconda::ucsc-liftover
```

- `bcftools`
- `GNU parallel`
- `Aldy4` (only required if running Aldy analysis)
- `UCSC LiftOver` (for reference genome conversions)

PGxIntegrator requires the following tools:

- **Python 3**
- `bcftools`
- `Aldy4`
- `Liftover` (if reference builds differ)
- `GNU parallel` (for batch processing)

## License

This project is licensed under the MIT License.

## Contributing

Pull requests and feature suggestions are welcome! Please create an issue for any bugs or feature requests.

## Acknowledgment

This project is developed as part of **Seoul National University, Graduate School of Data Science**.

## Contact

For questions or collaborations, contact **Chanhee Lee** at **achl1990@snu.ac.kr**.

---

🚀 **Automate your pharmacogenetic analysis with PGxIntegrator!**

