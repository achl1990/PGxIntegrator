# PGxIntegrator

**PGxIntegrator** is a toolkit for **integrating Whole Exome Sequencing (WES) and imputed Array data into a WGS-like VCF** for pharmacogenomic analysis. It also provides supplementary tools for **multi-sample VCF splitting, Liftover (hg19 → hg38), preprocessing, and parallel Aldy4-based star-allele genotyping**.

## Features
✅ **WES + Array Integration** → Merges WES and (imputed) array data into a single VCF.  
✅ **Multi-Sample VCF Processing** → Splits multi-sample VCFs into single samples and runs Aldy4 in parallel.  
✅ **Automated Preprocessing** → Standardizes, normalizes, and filters variants.  
✅ **Liftover Support** → Converts `hg19 → hg38` if needed.  
✅ **Flexible & Scalable** → Works with both single-sample and multi-sample VCFs.  

## Installation
PGxIntegrator requires Python and several bioinformatics tools, including `bcftools`, and `Liftover`. `Aldy4` is required only if running Aldy analysis.

```bash
# Clone the repository
git clone https://github.com/yourusername/PGxIntegrator.git
cd PGxIntegrator

# Install dependencies
pip install -r requirements.txt
```

## Usage

**⚠ Note:** The `--integrate` function is currently under review and will be disabled until the next update.
### **Basic Command**
```bash
python pgxintegrator.py --wes wes.vcf --array array.vcf --reference hs38DH.fa --gene CYP2C19
```

### **Full List of Arguments**

| Argument | Description | Default |
|----------|------------|---------|
| `--wes` | Path to WES VCF (or CRAM/BAM) | **Required** |
| `--array` | Path to imputed Array VCF | Optional |
| `--reference` | Path to reference genome (`hs38DH.fa`) | **Required** |
| `--reference-build-wes` | Build of WES data (`hg19` or `hg38`) | **Required** |
| `--reference-build-array` | Build of Array data (`hg19` or `hg38`) | **Required if --array** |
| `--liftover` | Apply hg19 → hg38 LiftOver | `False` |
| `--gene` | Gene(s) to analyze (`CYP2C19,VKORC1, all`) | **Required** |
| `--integrate` | Enable WES + Array integration (**Currently Disabled**) | `False` |
| `--harmonization` | Standardize, normalize, and filter | `True` |
| `--multi-sample` | Treat VCF as multi-sample | `False` |
| `--split-multi` | Split multi-sample VCF | `True` |
| `--run-aldy` | Run Aldy analysis | `True` |
| `--jobs` | Number of parallel jobs for GNU parallel | `8` |
| `--aldy-output` | Path to store Aldy results | `./results` |
| `--output-dir` | Directory for final files | `./results` |
| `--log` | Path to log file | `./log.txt` |
| `--keep-temp` | Retain intermediate files | `False` |
| `--jobs` | Number of CPU threads | `8` |
| `--dry-run` | Run pipeline without executing | `False` |

## Example Workflows
### **Integrate WES + Array and run Aldy**
```bash
python pgxintegrator.py --wes sample_wes.vcf --array sample_array.vcf \
  --reference hs38DH.fa --gene CYP2C19 --integrate --run-aldy
```

### **Split Multi-Sample VCF and Run Aldy in Parallel**
```bash
python pgxintegrator.py --wes multi_sample.vcf --split-multi --run-aldy --aldy-threads 16
```

### **Apply Liftover from hg19 to hg38**
```bash
python pgxintegrator.py --wes sample_wes_hg19.vcf --liftover --reference hs38DH.fa
```

## Dependencies
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

## Contact
For questions or collaborations, contact **[Your Name]** at **[Your Email]**.

---
🚀 **Automate your pharmacogenetic analysis with PGxIntegrator!**

