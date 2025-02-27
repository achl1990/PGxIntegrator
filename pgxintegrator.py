import argparse
import subprocess
import sys

parser = argparse.ArgumentParser(description="PGxIntegrator: Pharmacogenomic Analysis Pipeline")

# Required arguments
parser.add_argument("--wes", help="Path to WES VCF")
parser.add_argument("--array", help="Path to imputed Array VCF")
parser.add_argument("--reference", required=True, help="Reference genome (hs38DH.fa)")
parser.add_argument("--gene", required=True, help="Target genes for Aldy (comma-separated or 'all')")
parser.add_argument("--output-dir", required=True, help="Directory for output files")

# Optional arguments
parser.add_argument("--multi-sample", action="store_true", help="Indicates if input VCF is multi-sample")
parser.add_argument("--integrate", action="store_true", help="Enable WES + Array integration")    # Currently disabled
parser.add_argument("--run-aldy", action="store_true", help="Run Aldy4 for genotyping")
parser.add_argument("--process-aldy-solutions", action="store_true", help="Process Aldy4 solutions")
parser.add_argument("--visualize", action="store_true", help="Generate visualizations")

args = parser.parse_args()

# At least one of --wes or --array is required
if not args.wes and not args.array:
    print("Error: At least one of --wes or --array must be provided.", file=sys.stderr)
    sys.exit(1)

# Define output directories
extracted_vcf = f"{args.output_dir}/extracted_vcfs"
preprocessed_vcf = f"{args.output_dir}/preprocessed_vcfs"
integrated_vcf = f"{args.output_dir}/integrated_vcf"
single_sample_vcf = f"{args.output_dir}/single_sample_vcfs"
aldy_output = f"{args.output_dir}/aldy_results"
processed_aldy_results = f"{args.output_dir}/processed_aldy_results"
visualization_results = f"{args.output_dir}/visualization_results"

# Ensure output directory exists
subprocess.run(["mkdir", "-p", args.output_dir])

# Step 1: Extract VCF by Gene
input_vcf = args.wes if args.wes else args.array
subprocess.run(["bash", "extract_vcf_by_gene.sh", "--input", input_vcf, "--output", extracted_vcf])

# Step 2: Preprocess VCFs
subprocess.run(["bash", "preprocess_vcfs.sh", "--input", extracted_vcf, "--output", preprocessed_vcf])

# Step 3: Harmonization & Integration (Currently disabled for debugging)
if args.integrate:
    if not args.wes or not args.array:
        print("Error: Both --wes and --array are required for integration.", file=sys.stderr)
        sys.exit(1)
    subprocess.run(["python", "harmonize_integrate.py", "--wes", preprocessed_vcf, "--array", args.array, "--output", integrated_vcf])
    input_vcf = integrated_vcf  # Update input for next steps
else:
    input_vcf = preprocessed_vcf

# Step 4: Split Multi-Sample VCFs (if applicable)
if args.multi_sample:
    subprocess.run(["bash", "split_multi_sample_vcfs.sh", "--input", input_vcf, "--output", single_sample_vcf])
    input_vcf = single_sample_vcf  # Update input for next steps

# Step 5: Run Aldy in Parallel (if requested)
if args.run_aldy:
    subprocess.run(["bash", "run_aldy_parallel.sh", "--input", input_vcf, "--gene", args.gene, "--output", aldy_output])
    input_vcf = aldy_output  # Update input for next steps

# Step 6: Process Aldy Solutions (if requested)
if args.process_aldy_solutions:
    subprocess.run(["python", "process_aldy_solutions.py", "--input", aldy_output, "--output", processed_aldy_results])
    input_vcf = processed_aldy_results  # Update input for next steps

# Step 7: Visualization (if requested)
if args.visualize:
    subprocess.run(["python", "visualization.py", "--input", processed_aldy_results, "--output", visualization_results])

print("PGxIntegrator pipeline completed successfully.")
