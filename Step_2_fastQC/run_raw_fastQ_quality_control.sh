#!/usr/bin/env bash

#SBATCH --job-name=fastQ_files_quality_control_before_trimming
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=1000MB 
#SBATCH --time=2:00:00 
#SBATCH --partition=pibu_el8 
#SBATCH --array=1-30

# Output and error files go here:
# %A → main job ID
# %a → array task ID

#SBATCH --output=/data/users/apiatkowska/RNAseq/Step_2_fastQC/fastQ_results/output/output_%A_%a.o
#SBATCH --error=/data/users/apiatkowska/RNAseq/Step_2_fastQC/fastQ_results/error/error_%A_%a.e

# email notifications: send email when job starts and finishes

#SBATCH --mail-user=agnieszka.piatkowska@students.unibe.ch
#SBATCH --mail-type=begin,end

# Prevent automatic module load
module purge

# Set directories as variables:
# raw FASTQ input files
FASTQ_DIR=/data/courses/rnaseq_course/toxoplasma_de/reads_Blood
# saves quality control reports output files to: 
OUTPUT_DIR=/data/users/apiatkowska/RNAseq/Step_2_fastQC/fastQC_results
# path to FastQC apptainer:
CONTAINER=/containers/apptainer/fastqc-0.12.1.sif

# Creates output directory if it doesn’t exist
# -p prevents error if it already exists
mkdir -p "$OUTPUT_DIR"

# collect fastq files
# Only select *.fastq.gz files, ignore previous FastQC outputs
FASTQ_FILES=($(find "$FASTQ_DIR" -maxdepth 1 -type f -name "*.fastq.gz" ! -name "*_fastqc*"))
# counts how many fastq files exists
NUM_FILES=${#FASTQ_FILES[@]}

# select the fastq file for this array task
# checks “Is this array task ID valid for the number of files?”
# Prevents crashes if array size > number of files

if [ $SLURM_ARRAY_TASK_ID -le $NUM_FILES ]; then
    FASTQ_FILE=${FASTQ_FILES[$SLURM_ARRAY_TASK_ID-1]}
    # Prints message to output file
    echo "Running FastQC on $FASTQ_FILE"
    # Run FastQC inside the container
    # --bind mounts a directory from the host into the container
    apptainer exec --bind $FASTQ_DIR:$FASTQ_DIR --bind $OUTPUT_DIR:$OUTPUT_DIR "$CONTAINER" fastqc --outdir "$OUTPUT_DIR" "$FASTQ_FILE"
# If task ID > number of FASTQ files
# job does nothing instead of crashing
else
    echo "Task $SLURM_ARRAY_TASK_ID skipped: no corresponding FASTQ file"
fi