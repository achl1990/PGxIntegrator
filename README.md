# PGxIntegrator

**PGxIntegrator** is a **pipeline** for **integrating Whole Exome Sequencing (WES) and imputed Array data into a WGS-like VCF** for pharmacogenomic analysis. It is a pipelining tool that integrates various bioinformatics tools for ease of use, offering **multi-sample VCF splitting, Liftover (hg19 → hg38), preprocessing, and parallel Aldy4-based star-allele genotyping**.


## How PGxIntegrator Works

PGxIntegrator **integrates Whole Exome Sequencing (WES) and imputed Array data** to create a **WGS-like dataset** for pharmacogenomic analysis. This enables better variant calling and star allele genotyping.

The diagram below illustrates the integration process:

![WES + Array Integration](https://raw.githubusercontent.com/achl1990/PGxIntegrator/main/wes_array_integration_figure.png)

- **WES (Top Line):** Captures exonic regions but misses intronic or regulatory variants.
- **Array (Middle Line):** Captures known pharmacogenomic SNPs and structural variants.
- **WES + Array (Bottom Line):** Combines both sources, ensuring comprehensive variant detection.

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
| `--wes` | Path to WES VCF (or CRAM/BAM) | **At least one of --wes or --array is required. You may provide both.** |
| `--array` | Path to imputed Array VCF | **At least one of --wes or --array is required. You may provide both.** |
| `--fetch-coordinates` | Fetch gene coordinates from Ensembl and save to a file | `False` |
| `--build` | Genome build for Ensembl coordinate fetching (`hg19` or `hg38`) | `hg38` |
| `--buffer-size` | Buffer size around gene coordinates (in bp) | `10000` |
| `--reference`             | Path to reference genome (`hs38DH.fa`)                  | **Required**            |
| `--gene`                  | Gene(s) to analyze (`CYP2C19,VKORC1, all`)              | **Required**            |
| `--output-dir`            | Directory for final files                               | **Required**            |
| `--integrate`             | Enable WES + Array integration (**Currently Disabled**) | `False`                 |
| `--multi-sample`          | Treat VCF as multi-sample                               | `False`                 |
| `--split-multi`           | Split multi-sample VCF                                  | `True`                  |
| `--run-aldy`              | Run Aldy analysis                                       | `True`                  |
| `--jobs`                  | Number of parallel jobs for GNU parallel                | `8`                     |
| `--aldy-output`           | Path to store Aldy results                              | `./results`             |
| `--log`                   | Path to log file                                        | `./log.txt`             |
| `--keep-temp`             | Retain intermediate files                               | `False`                 |
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

### **Fetch Updated Coordinates for Pharmacogenes**
```bash
python pgxintegrator.py --fetch-coordinates --build hg38 --buffer-size 10000 --output-dir results
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

