#!/usr/bin/env bash
#SBATCH --job-name=multiqc_fastp
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/apiatkowska/RNAseq/Step_2.2_trimmed_fastQ/multiQC_fastp_results/multiqc_%j.out
#SBATCH --error=/data/users/apiatkowska/RNAseq/Step_2.2_trimmed_fastQ/multiQC_fastp_results/multiqc_%j.err
#SBATCH --mail-user=agnieszka.piatkowska@students.unibe.ch
#SBATCH --mail-type=BEGIN,END

# Set directories as variables:
# trimmed_fastQ files serve as input files
INPUT_DIR=/data/users/apiatkowska/RNAseq/Step_2.2_trimmed_fastQ/trimmed_fastQ
# for multiQC report
OUTPUT_DIR=/data/users/apiatkowska/RNAseq/Step_2.2_trimmed_fastQ/multiQC_fastp_results
# for container
CONTAINER=/containers/apptainer/multiqc-1.19.sif

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Bind trimmed FASTQ folder to a simple path inside container
apptainer exec --bind "$INPUT_DIR:/data_trimmed_fastq" \
    "$CONTAINER" \
    multiqc /data_trimmed_fastq \
    --outdir "$OUTPUT_DIR" \
    --filename "multiqc_trimmed_fastq_report.html" \
    --force