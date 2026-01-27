#!/usr/bin/env bash

#SBATCH --job-name=multiqc_on_raw_fastQ_files
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --partition=pibu_el8
#SBATCH --output=multiqc.out
#SBATCH --error=multiqc.err

#SBATCH --output=/data/users/apiatkowska/RNAseq/Step_2_fastQC/multiQC_results/multiqc_%j.out
#SBATCH --error=/data/users/apiatkowska/RNAseq/Step_2_fastQC/multiQC_results/multiqc_%j.err

#SBATCH --mail-user=agnieszka.piatkowska@students.unibe.ch
#SBATCH --mail-type=BEGIN,END

# Set directories as variables:
# takes input from
FASTQC_DIR=/data/users/apiatkowska/RNAseq/Step_2_fastQC/fastQC_results
# writes output here
MULTIQC_DIR=/data/users/apiatkowska/RNAseq/Step_2_fastQC/multiQC_results
# (path to MultiQC apptainer) takes aptteiner from:
CONTAINER=/containers/apptainer/multiqc-1.19.sif

# Make sure output directory exists
mkdir -p "$MULTIQC_DIR"

# Run MultiQC
apptainer exec \
  --bind "$FASTQC_DIR:$FASTQC_DIR" \
  --bind "$MULTIQC_DIR:$MULTIQC_DIR" \
  "$CONTAINER" \
  multiqc "$FASTQC_DIR" --outdir "$MULTIQC_DIR"