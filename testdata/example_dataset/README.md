Example dataset for PGxIntegrator
================================

Structure
---------
wes/
  - chr10_wes_hg38.vcf.gz
  - VKORC1_wes_hg38.vcf.gz
array/
  - chr10_array_hg19.vcf.gz
  - VKORC1_array_hg19.vcf.gz
manifests/
  - sample_manifest.tsv
  - gene_list.tsv
  - dataset_manifest.tsv
configs/
  - example_config.yaml
outputs/
  - prepared/
  - split/
  - aldy/

Design
------
- Mixed layout by design:
  * CYP2C19 uses chromosome-format inputs
  * VKORC1 uses gene-format inputs
- 5 samples with different WES and array sample IDs
- WES is hg38 and includes exonic variants only
- Array is hg19 and includes important non-exonic plus some exonic variants
- Synthetic hg19 array REF alleles are made liftOver-consistent with hg38
  so bcftools norm succeeds after liftOver.
