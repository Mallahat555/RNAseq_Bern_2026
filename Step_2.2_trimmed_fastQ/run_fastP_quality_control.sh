#!/usr/bin/env bash

#SBATCH --job-name=individual_fastP_quality_control
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=1000MB 
#SBATCH --time=1:00:00 
#SBATCH --partition=pibu_el8 
#SBATCH --array=1-30

#SBATCH --output=/data/users/apiatkowska/RNAseq/Step_2.2_trimmed_fastQ/individual_fastP_quality_control_reports/output/output_%A_%a.o
#SBATCH --error=/data/users/apiatkowska/RNAseq/Step_2.2_trimmed_fastQ/individual_fastP_quality_control_reports/error/error_%A_%a.e

#SBATCH --mail-user=agnieszka.piatkowska@students.unibe.ch
#SBATCH --mail-type=BEGIN,END

# Prevent automatic module load
module purge

# Set directories as variables:
# takes input from: 
INPUT_DIR=/data/users/apiatkowska/RNAseq/Step_2.2_trimmed_fastQ/trimmed_fastQ
# saves output to: 
OUTPUT_DIR=/data/users/apiatkowska/RNAseq/Step_2.2_trimmed_fastQ/individual_fastP_quality_control_reports
# path to FastQC apptainer
CONTAINER=/containers/apptainer/fastqc-0.12.1.sif

# Creates output directory if it doesn’t exist
# -p prevents error if it already exists
mkdir -p "$OUTPUT_DIR"/output
mkdir -p "$OUTPUT_DIR"/error

# Collect ONLY fastp-trimmed FASTQ files
FASTQ_FILES=($(ls "$INPUT_DIR"/*_trimmed.fastq.gz))
NUM_FILES=${#FASTQ_FILES[@]}

# Select file for this array task
if [ "$SLURM_ARRAY_TASK_ID" -le "$NUM_FILES" ]; then
    FASTQ_FILE=${FASTQ_FILES[$SLURM_ARRAY_TASK_ID-1]}
    echo "Running FastQC on fastp-trimmed file: $FASTQ_FILE"

    apptainer exec \
        --bind /data \
        "$CONTAINER" \
        fastqc \
        --outdir "$OUTPUT_DIR" \
        "$FASTQ_FILE"
else
    echo "No FASTQ file for array task $SLURM_ARRAY_TASK_ID"
fi


