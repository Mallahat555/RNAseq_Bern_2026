#!/bin/bash
#SBATCH --job-name=extract_igv_region
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --partition=pibu_el8
#SBATCH --array=0-14

#SBATCH --output=/data/users/apiatkowska/RNAseq/Step_3_map_reads/mapped_reads_small_regions_IGV/logs/%x_%A_%a.out
#SBATCH --error=/data/users/apiatkowska/RNAseq/Step_3_map_reads/mapped_reads_small_regions_IGV/logs/%x_%A_%a.err

#SBATCH --mail-user=agnieszka.piatkowska@students.unibe.ch
#SBATCH --mail-type=BEGIN,END

# Samples list
SAMPLES=(
  SRR7821949 SRR7821950 SRR7821951 SRR7821952 SRR7821953
  SRR7821954 SRR7821955 SRR7821956 SRR7821957 SRR7821968
  SRR7821969 SRR7821970 SRR7821971 SRR7821972 SRR7821973
)

# one sample is processed by one array
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Extract small region for IGV
REGION="chr1:150000000-150100000"
CHR=$(echo "$REGION" | cut -d':' -f1)
START=$(echo "$REGION" | cut -d':' -f2 | cut -d'-' -f1)
END=$(echo "$REGION" | cut -d':' -f2 | cut -d'-' -f2)
REGION_TAG="${CHR}_${START}_${END}"
echo "Processing $SAMPLE | Region: $REGION"

# Set directories as variables:
# Apptainer container with samtools
CONTAINER=/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif
# Directory containing input coordinate-sorted BAM files
INPUT_BAM_DIR=/data/users/apiatkowska/RNAseq/Step_3_map_reads/mapped_reads/bam_sorted
# Base output directory for BAM files restricted to small genomic regions (IGV)
OUTPUT_BASE=/data/users/apiatkowska/RNAseq/Step_3_map_reads/mapped_reads_small_regions_IGV
# Directory for region-specific sorted BAM output files
OUTPUT_BAM_DIR=$OUTPUT_BASE/bam_sorted_small_reads
# Directory for standard output and error log files
LOG_DIR=$OUTPUT_BASE/logs
# Input BAM file for the current sample (coordinate-sorted)
INPUT_BAM=${INPUT_BAM_DIR}/${SAMPLE}_sorted.bam
# Output BAM file containing reads from the specified genomic region
OUTPUT_BAM=${OUTPUT_BAM_DIR}/${SAMPLE}_${REGION_TAG}.bam

# Create directories
mkdir -p "$OUTPUT_BAM_DIR"
mkdir -p "$LOG_DIR"

# Extract region
apptainer exec --bind /data "$CONTAINER" \
  samtools view -b "$INPUT_BAM" "$REGION" > "$OUTPUT_BAM"

# Index BAM
apptainer exec --bind /data "$CONTAINER" \
  samtools index "$OUTPUT_BAM"

echo "Finished $SAMPLE ($REGION_TAG)"

