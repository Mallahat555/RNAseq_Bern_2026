#!/usr/bin/env bash
#SBATCH --job-name=trimming_fastQ_files_pair_ended
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8
#SBATCH --array=1-30

#SBATCH --output=/data/users/apiatkowska/RNAseq/Step_2.2_trimmed_fastQ/trimmed_fastQ/output/output_%A_%a.o
#SBATCH --error=/data/users/apiatkowska/RNAseq/Step_2.2_trimmed_fastQ/trimmed_fastQ/error/error_%A_%a.e

#SBATCH --mail-user=agnieszka.piatkowska@students.unibe.ch
#SBATCH --mail-type=BEGIN,END

# Set directories as variables:
FASTQ_DIR=/data/courses/rnaseq_course/toxoplasma_de/reads_Blood      # raw FASTQ input files
OUTPUT_DIR=/data/users/apiatkowska/RNAseq/Step_2.2_trimmed_fastQ/trimmed_fastQ  # folder for trimmed fastp output files
CONTAINER=/containers/apptainer/fastp_0.24.1.sif   # fastp tool from container

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# List of R1 files (forward reads)
FASTQ_FILES_R1=($(ls ${FASTQ_DIR}/*_1.fastq.gz | grep -v '_trimmed'))
NUM_FILES=${#FASTQ_FILES_R1[@]}

# Pick file for this array task 
if [ $SLURM_ARRAY_TASK_ID -le $NUM_FILES ]; then
    FASTQ_FILE_R1=${FASTQ_FILES_R1[$SLURM_ARRAY_TASK_ID-1]}
    BASENAME=$(basename "$FASTQ_FILE_R1" _1.fastq.gz)
    FASTQ_FILE_R2="${FASTQ_DIR}/${BASENAME}_2.fastq.gz"

    # Output trimmed files
    TRIMMED_R1="${OUTPUT_DIR}/${BASENAME}_1_trimmed.fastq.gz"
    TRIMMED_R2="${OUTPUT_DIR}/${BASENAME}_2_trimmed.fastq.gz"

    echo "Trimming paired-end reads:"
    echo "R1: $FASTQ_FILE_R1 -> $TRIMMED_R1"
    echo "R2: $FASTQ_FILE_R2 -> $TRIMMED_R2"

# Run fastp in paired-end mode
    apptainer exec \
      --bind "$FASTQ_DIR:$FASTQ_DIR" \
      --bind "$OUTPUT_DIR:$OUTPUT_DIR" \
      "$CONTAINER" \
      fastp \
        -i "$FASTQ_FILE_R1" \
        -I "$FASTQ_FILE_R2" \
        -o "$TRIMMED_R1" \
        -O "$TRIMMED_R2" \
        -h "$OUTPUT_DIR/${BASENAME}_fastp.html" \
        -j "$OUTPUT_DIR/${BASENAME}_fastp.json" \
        -w $SLURM_CPUS_PER_TASK
else
    echo "Task $SLURM_ARRAY_TASK_ID skipped: no corresponding FASTQ pair"
fi