import subprocess
import os
import sys
import logging
import random
import time

# Generate a unique identifier based on timestamp and random number
unique_id = f"{int(time.time())}_{random.randint(1000, 9999)}"

# Setup logging with unique ID to avoid overwriting
log_file = f"combined_variants_{unique_id}.log"
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger()

# Check for the correct number of arguments
if len(sys.argv) != 3:
    logger.error("Usage: python combined_variants.py <working_directory> <sample_info_file>")
    sys.exit(1)

# Command-line arguments
project_dir = sys.argv[1]
sample_info_file = sys.argv[2]

# Define directories
job_files_dir = os.path.join(project_dir, "job_files")
job_output_dir = os.path.join(project_dir, "job_output")
job_error_dir = os.path.join(project_dir, "job_error")
mutect2_dir = os.path.join(project_dir, "mutect2")
freebayes_dir = os.path.join(project_dir, "freebayes")
merged_variants_dir = os.path.join(project_dir, "merged_variants")

# Ensure directories exist
os.makedirs(job_files_dir, exist_ok=True)
os.makedirs(job_output_dir, exist_ok=True)
os.makedirs(job_error_dir, exist_ok=True)
os.makedirs(mutect2_dir, exist_ok=True)
os.makedirs(freebayes_dir, exist_ok=True)
os.makedirs(merged_variants_dir, exist_ok=True)

# SLURM job templates
mutect2_job_template = """#!/bin/bash
#SBATCH --job-name=mutect2_{sample}.%j
#SBATCH --output={output_dir}/{sample}_mutect2.out.%j
#SBATCH --error={error_dir}/{sample}_mutect2.err.%j
#SBATCH --ntasks=2
#SBATCH --nodes=1

docker run --rm -i \
    -v {project_dir}:/data \
    -v /path/to/ref/hg38:/ref \
    -v /database/mutect_ref:/panel_dir \
    -v /database/mutect_ref:/germline_dir \
    /{docker_container} /bin/bash -c "
        mkdir -p /data/mutect2 && chmod 777 /data/mutect2
        echo \$(date)': Starting samtools addreplacerg.' | tee /dev/stderr

        samtools addreplacerg --threads 6 -r ID:{sample} -r SM:{sample} -o /data/bam/{sample}.filtered.recalibrated.rg.bam /data/bam/{{sample}.umi-processed_dupes-removed.bam

        samtools index --threads 6 /data/bam/{sample}.filtered.recalibrated.rg.bam

        /root/gatk-4.3.0.0/gatk Mutect2 -R /ref/hg38.fa -I /data/bam/{sample}.filtered.recalibrated.rg.bam -L /ref/{interval_list} --tumor-lod-to-emit 0.0 --minimum-allele-fraction 0.0 --f1r2-tar-gz /data/mutect2/{sample}_f1r2.tar.gz --panel-of-normals /panel_dir/1000g_pon.hg38.vcf.gz --germline-resource /germline_dir/af-only-gnomad.hg38.vcf.gz -O /data/mutect2/variants_{sample}.vcf.gz
        /root/gatk-4.3.0.0/gatk FilterMutectCalls -R /ref/hg38.fa -V /data/mutect2/variants_{sample}.vcf.gz --orientation-bias-artifact-priors /data/mutect2/{sample}_artifact-prior.tar.gz -O /data/mutect2/variants_{sample}.filtered.vcf.gz
        rm /data/bam/{sample}.filtered.recalibrated.rg.bam
    "
"""

freebayes_job_template = """#!/bin/bash
#SBATCH --job-name=freebayes_{sample}.%j
#SBATCH --output={output_dir}/{sample}_freebayes.out.%j
#SBATCH --error={error_dir}/{sample}_freebayes.err.%j
#SBATCH --ntasks=2
#SBATCH --nodes=1

{freebayes_path} --fasta-reference hg38.fa --bam {project_dir}/bam/{sample}.umi-processed_dupes-removed.bam \
    --targets {bed_file} --min-alternate-fraction 0.0001 --min-alternate-count 5 --pooled-discrete --ploidy 2 \
    > {freebayes_dir}/freebayes_{sample}.vcf

# Sorting, compressing, indexing, and merging
{bcftools_path} sort -o {freebayes_dir}/{sample}.freebayes.sorted.vcf {freebayes_dir}/freebayes_{sample}.vcf
{bgzip_path} {freebayes_dir}/{sample}.freebayes.sorted.vcf
{tabix_path} {freebayes_dir}/{sample}.freebayes.sorted.vcf.gz

# Merging with Mutect2 VCF
{bcftools_path} concat -a {mutect2_dir}/variants_{sample}.filtered.vcf.gz {freebayes_dir}/{sample}.freebayes.sorted.vcf.gz \
    -o {merged_variants_dir}/merged_{sample}.vcf.gz
{bcftools_path} norm -d none {merged_variants_dir}/merged_{sample}.vcf.gz -o {merged_variants_dir}/normalized_{sample}.vcf
"""

# VEP job template
vep_job_template = """#!/bin/bash
#SBATCH --job-name=vep_{sample}.%j
#SBATCH --output={output_dir}/{sample}_vep.out.%j
#SBATCH --error={error_dir}/{sample}_vep.err.%j
#SBATCH --ntasks=2
#SBATCH --nodes=1

sudo chmod -R 777 /data/merged_variants

docker run --rm -i \
    -v {project_dir}:/data \
    -v /path/to/ref/hg38:/ref \
    -v /path/to/vep_data:/vep_data \
    -v /path/to/VEP_plugins:/plugins ensemblorg/ensembl-vep \
    /bin/bash -c "
        echo \$(date)': Starting VEP annotation.' | tee /dev/stderr
        vep --input_file /data/merged_variants/normalized_{sample}.vcf --output_file /data/merged_variants/variants_{sample}.vep.vcf --format vcf --vcf --symbol --terms SO --tsl --biotype \
        --hgvs --fasta /ref/hg38.fa --offline --cache --dir_cache /vep_data --everything --dir_plugins /plugins --force_overwrite 
    "
"""

# Funcotator job template
funcotator_job_template = """#!/bin/bash
#SBATCH --job-name=funcotator_{sample}.%j
#SBATCH --output={output_dir}/{sample}_funcotator.out.%j
#SBATCH --error={error_dir}/{sample}_funcotator.err.%j
#SBATCH --ntasks=2
#SBATCH --nodes=1

docker run --rm -i \
    -v {project_dir}:/data \
    -v /path/to/ref/hg38:/ref \
    -v /path/to/funcotator_dataSources.v1.7.20200521s:/func_data_source \
    {docker_container} /bin/bash -c "
        echo \$(date)': Starting Funcotator.' | tee /dev/stderr
        /root/gatk-4.3.0.0/gatk Funcotator --variant /data/merged_variants/normalized_{sample}.vcf --reference /ref/hg38.fa --ref-version hg38 \
        --data-sources-path /func_data_source --output /data/merged_variants/variants_{sample}.funcotated.vcf --output-file-format VCF --disable-sequence-dictionary-validation 
        echo \$(date)': Funcotator complete.'
    "
"""

# Paths for FreeBayes and BCFtools
freebayes_path = "/path/to/freebayes"
bcftools_path = "/path/to/bcftools"
bgzip_path = "/path/to/bgzip"
tabix_path = "/path/to/tabix"

# Store Mutect2 job IDs
mutect2_job_ids = []

# Read the sample info file and generate job files
with open(sample_info_file, 'r') as f:
    next(f)  # Skip header line
    for line in f:
        sample = line.strip().split('\t')[0]  # Assuming the first column is sample ID

        # Generate and write Mutect2 job file
        mutect2_job_file_content = mutect2_job_template.format(
            sample=sample,
            output_dir=job_output_dir,
            error_dir=job_error_dir,
            project_dir=project_dir
        )
        mutect2_job_file_path = os.path.join(job_files_dir, f"{sample}_mutect2.sh")
        with open(mutect2_job_file_path, 'w') as job_file:
            job_file.write(mutect2_job_file_content)

        # Submit Mutect2 job and store the job ID
        mutect2_cmd = f"sbatch --parsable {mutect2_job_file_path}"
        logger.info(f"Submitting Mutect2 job for {sample}")
        result = subprocess.run(mutect2_cmd, shell=True, capture_output=True, text=True)

        if result.returncode != 0:
            logger.error(f"Error submitting Mutect2 job for {sample}: {result.stderr}")
        else:
            job_id = result.stdout.strip()
            mutect2_job_ids.append(job_id)
            logger.info(f"Mutect2 job submitted for {sample} with Job ID: {job_id}")
            logger.info(f"Mutect2 job running for {sample}: Job ID {job_id}")

# Create a dependency string for FreeBayes jobs
mutect2_dependency = ":".join(mutect2_job_ids)

# Store freebayes job IDs
freebayes_job_ids = []

# Submit FreeBayes jobs with dependency on Mutect2 completion
with open(sample_info_file, 'r') as f:
    next(f)  # Skip header line
    for line in f:
        sample = line.strip().split('\t')[0]

        # Generate and write FreeBayes job file
        freebayes_job_file_content = freebayes_job_template.format(
            sample=sample,
            output_dir=job_output_dir,
            error_dir=job_error_dir,
            project_dir=project_dir,
            freebayes_dir=freebayes_dir,
            merged_variants_dir=merged_variants_dir,
            freebayes_path=freebayes_path,
            bcftools_path=bcftools_path,
            bgzip_path=bgzip_path,
            tabix_path=tabix_path,
            mutect2_dir=mutect2_dir
        )
        freebayes_job_file_path = os.path.join(job_files_dir, f"{sample}_freebayes.sh")
        with open(freebayes_job_file_path, 'w') as job_file:
            job_file.write(freebayes_job_file_content)

        # Submit FreeBayes job with dependency on Mutect2 jobs
        freebayes_cmd = f"sbatch --dependency=afterok:{mutect2_dependency} --parsable {freebayes_job_file_path}"
        logger.info(f"Submitting FreeBayes job for {sample}")
        freebayes_result = subprocess.run(freebayes_cmd, shell=True, capture_output=True, text=True)

        if freebayes_result.returncode != 0:
            logger.error(f"Error submitting FreeBayes job for {sample}: {freebayes_result.stderr}")
        else:
            job_id = freebayes_result.stdout.strip()
            freebayes_job_ids.append(job_id)
            logger.info(f"FreeBayes job submitted for {sample} with Job ID: {job_id}")
            logger.info(f"FreeBayes job running for {sample}: Job ID {job_id}")

# VEP and Funcotator submission
vep_job_ids = []
funcotator_job_ids = []

# Create a dependency string for annotation jobs
variant_calling_dependency = ":".join(freebayes_job_ids)

with open(sample_info_file, 'r') as f:
    next(f)  # Skip header line
    for line in f:
        sample = line.strip().split('\t')[0]  # Assuming the first column is sample ID

        # Generate and write VEP job file
        vep_job_file_content = vep_job_template.format(
            sample=sample,
            output_dir=job_output_dir,
            error_dir=job_error_dir,
            project_dir=project_dir
        )
        vep_job_file_path = os.path.join(job_files_dir, f"{sample}_vep.sh")
        with open(vep_job_file_path, 'w') as job_file:
            job_file.write(vep_job_file_content)

        # Submit VEP job with dependency on FreeBayes jobs
        vep_cmd = f"sbatch --dependency=afterok:{variant_calling_dependency} --parsable {vep_job_file_path}"
        vep_result = subprocess.run(vep_cmd, shell=True, capture_output=True, text=True)

        if vep_result.returncode != 0:
            logger.error(f"Error submitting VEP job for {sample}: {vep_result.stderr}")
        else:
            vep_job_id = vep_result.stdout.strip()
            vep_job_ids.append(vep_job_id)
            logger.info(f"VEP job submitted for {sample} with Job ID: {vep_job_id}")

        # Generate and write Funcotator job file
        funcotator_job_file_content = funcotator_job_template.format(
            sample=sample,
            output_dir=job_output_dir,
            error_dir=job_error_dir,
            project_dir=project_dir
        )
        funcotator_job_file_path = os.path.join(job_files_dir, f"{sample}_funcotator.sh")
        with open(funcotator_job_file_path, 'w') as job_file:
            job_file.write(funcotator_job_file_content)

        # Submit Funcotator job with dependency on VEP jobs
        funcotator_cmd = f"sbatch --dependency=afterok:{vep_job_id} --parsable {funcotator_job_file_path}"
        funcotator_result = subprocess.run(funcotator_cmd, shell=True, capture_output=True, text=True)

        if funcotator_result.returncode != 0:
            logger.error(f"Error submitting Funcotator job for {sample}: {funcotator_result.stderr}")
        else:
            funcotator_job_id = funcotator_result.stdout.strip()
            funcotator_job_ids.append(funcotator_job_id)
            logger.info(f"Funcotator job submitted for {sample} with Job ID: {funcotator_job_id}")

logger.info("All jobs have been submitted.")

# Function to check job status
def check_job_status(job_id):
    job_status_cmd = f"scontrol show job {job_id} | grep JobState"
    result = subprocess.run(job_status_cmd, shell=True, capture_output=True, text=True)
    if "COMPLETED" in result.stdout:
        logger.info(f"Job {job_id} completed successfully.")
    elif "FAILED" in result.stdout:
        logger.error(f"Job {job_id} failed.")
    else:
        logger.info(f"Job {job_id} is still running or pending.")
        
# Example usage: call check_job_status with each job ID as needed.
for job_id in mutect2_job_ids:
    check_job_status(job_id)
