# ctDNA Variant Calling and Annotation Pipeline

This repository contains a Python-based workflow for **variant calling and annotation of ctDNA samples** using **SLURM**, **Docker**, and multiple industry-standard tools:

- **GATK Mutect2** – somatic variant calling  
- **FreeBayes** – additional low-frequency variant calling  
- **bcftools** – VCF processing and merging  
- **VEP** – variant annotation  
- **GATK Funcotator** – functional and clinical annotation  

The script automates job creation, submission, dependency handling, and logging for large-scale ctDNA analyses on an HPC cluster.

---

## Overview of the Pipeline

For each sample, the pipeline performs the following steps:

1. **Mutect2**
   - Adds read groups
   - Calls somatic variants
   - Filters variants using GATK best practices

2. **FreeBayes**
   - Calls variants with low allele fractions
   - Sorts, compresses, and indexes VCFs

3. **Variant Merging**
   - Merges Mutect2 and FreeBayes calls
   - Normalizes merged VCFs

4. **Annotation**
   - **VEP**: Comprehensive variant annotation
   - **Funcotator**: Functional and clinical annotation

Each step is submitted as a separate SLURM job with appropriate job dependencies.

---

## Requirements

### System
- Linux HPC environment
- SLURM workload manager
- Docker
- Python ≥ 3.7

### Software / Tools
- GATK 4.x (Dockerized)
- samtools
- FreeBayes
- bcftools
- bgzip
- tabix
- Ensembl VEP (Docker image)
- Reference genome: **hg38**

### Databases
- Mutect2 Panel of Normals (PoN)
- gnomAD germline resource
- VEP cache (offline)
- Funcotator data sources

---

## Directory Structure

The script creates and uses the following directory structure inside the project directory:

```
project_dir/
├── bam/                     # Input BAM files
├── job_files/               # Generated SLURM job scripts
├── job_output/              # SLURM stdout files
├── job_error/               # SLURM stderr files
├── mutect2/                 # Mutect2 outputs
├── freebayes/               # FreeBayes outputs
├── merged_variants/         # Merged and annotated VCFs
├── combined_variants.log    # Pipeline log file
```

---

## Input Files

### 1. Sample Info File

A tab-delimited text file with a header.  
The **first column must contain the sample ID**.

Example:

```
sample_id	other_metadata
SAMPLE_001	...
SAMPLE_002	...
```

---

### 2. BAM Files

Expected BAM locations:

- **Mutect2 input**
- **FreeBayes input**
  ```
  bam/{sample}.umi-processed_dupes-removed.bam
  ```

---

## Configuration

Before running the script, update the following paths in the Python file:

```python
freebayes_path = "/path/to/freebayes"
bcftools_path  = "/path/to/bcftools"
bgzip_path     = "/path/to/bgzip"
tabix_path     = "/path/to/tabix"
```

Also update Docker mount paths for:

- hg38 reference genome  
- Mutect2 Panel of Normals and germline resources  
- VEP cache and plugins  
- Funcotator data sources  

---

## Usage

Run the script from the command line:

```bash
python variant_Cal.py <project_dir> <sample_info_file>
```

### Example

```bash
python variant_calling_and_annotation_ctDNA.py /data/ctdna_project samples.tsv
```

---

## SLURM Job Dependencies

The workflow enforces the following execution order:

```
Mutect2
   ↓
FreeBayes + Merge
   ↓
VEP
   ↓
Funcotator
```

Dependencies are managed using:

```bash
sbatch --dependency=afterok
```

---

## Logging

A **unique log file** is created for each run:

```
combined_variants_<timestamp>_<random>.log
```

Logs include:
- Job submission status
- SLURM job IDs
- Error messages
- Job completion checks

---

## Job Monitoring

The script includes a helper function to check job status:

```python
check_job_status(job_id)
```

This uses `scontrol show job` to report:

- `COMPLETED`
- `FAILED`
- `RUNNING` / `PENDING`

---

## Outputs

Final annotated files are written to:

```
merged_variants/
├── normalized_<sample>.vcf
├── variants_<sample>.vep.vcf
├── variants_<sample>.funcotated.vcf
```

---

## Notes & Best Practices

- Ensure all reference files match **hg38**
- Confirm Docker images are accessible on compute nodes
- Validate SLURM resource requests (`--ntasks`, `--nodes`) for your cluster
- Test with a single sample before large batch runs

---

## Disclaimer

This pipeline is intended for **research use only**.  
Clinical applications require validation, quality control, and regulatory compliance.

---

## Author

Developed for automated ctDNA variant calling and annotation in HPC environments.
