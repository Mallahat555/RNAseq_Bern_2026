#!/bin/bash
#SBATCH --job-name=mapping_hisat2_for_all_samples
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --time=1:00:00
#SBATCH --partition=pibu_el8
#SBATCH --array=0-14

#SBATCH --output=/data/users/apiatkowska/RNAseq/Step_3_map_reads/mapped_reads/logs/%x_%A_%a.out
#SBATCH --error=/data/users/apiatkowska/RNAseq/Step_3_map_reads/mapped_reads/logs/%x_%A_%a.err

#SBATCH --mail-user=agnieszka.piatkowska@students.unibe.ch
#SBATCH --mail-type=BEGIN,END

# Sample list
SAMPLES=(
  SRR7821949
  SRR7821950
  SRR7821951
  SRR7821952
  SRR7821953
  SRR7821954
  SRR7821955
  SRR7821956
  SRR7821957
  SRR7821968
  SRR7821969
  SRR7821970
  SRR7821971
  SRR7821972
  SRR7821973
)

# one sample is processed by one array
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "Processing sample: $SAMPLE"

# set directories as variables
# Apptainer container with HISAT2 and samtools
CONTAINER=/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif
# Directory containing HISAT2 index files for the mouse reference genome (mm39)
INDEX_DIR=/data/users/apiatkowska/RNAseq/Step_3_map_reads/hisat2_indexing/index/mm39
# Directory containing trimmed paired-end FASTQ files
FASTQ_DIR=/data/users/apiatkowska/RNAseq/Step_2.2_trimmed_fastQ/trimmed_fastQ
# Base output directory for read mapping results
OUT_DIR=/data/users/apiatkowska/RNAseq/Step_3_map_reads/mapped_reads
# Directory for SAM alignment files
SAM_DIR=$OUT_DIR/sam
# Directory for unsorted BAM files
BAM_DIR=$OUT_DIR/bam
# Directory for coordinate-sorted BAM files
SORTED_DIR=$OUT_DIR/bam_sorted
# Directory for standard output and error log files
LOG_DIR=$OUT_DIR/logs

# Create directory structure
mkdir -p $OUT_DIR/sam
mkdir -p $OUT_DIR/bam
mkdir -p $OUT_DIR/bam_sorted
mkdir -p $OUT_DIR/logs

# Input FASTQ files
FASTQ_R1=${FASTQ_DIR}/${SAMPLE}_1_trimmed.fastq.gz
FASTQ_R2=${FASTQ_DIR}/${SAMPLE}_2_trimmed.fastq.gz


# Read mapping with HISAT2 (strand-specific RF)
# Print status message
echo "Running HISAT2 for $SAMPLE"

# Execute HISAT2 inside the Apptainer container
apptainer exec --bind /data "$CONTAINER" hisat2 \
  -x "$INDEX_DIR" \                  # HISAT2 index for the reference genome (mm39)
  -1 "$FASTQ_R1" \                   # Forward reads (R1) FASTQ file
  -2 "$FASTQ_R2" \                   # Reverse reads (R2) FASTQ file
  --rna-strandness RF \              # Strand-specific RNA-seq library (reverse-forward)
  -p "$SLURM_CPUS_PER_TASK" \        # Number of CPU threads to use
  -S "$SAM_DIR/${SAMPLE}.sam"        # Output SAM alignment file


# convert SAM to BAM
# Print status message
echo "Converting SAM to BAM for $SAMPLE"

# Convert SAM alignment file to BAM format using samtools
apptainer exec --bind /data "$CONTAINER" \
  samtools view \
    -@ "$SLURM_CPUS_PER_TASK" \       # Number of CPU threads
    -bS \                             # Output BAM format; input is SAM
    "$SAM_DIR/${SAMPLE}.sam" \        # Input SAM file
    > "$BAM_DIR/${SAMPLE}.bam"        # Output BAM file


# Sort BAM by genomic coordinates
# Print status message
echo "Sorting BAM for $SAMPLE"

# Sort BAM file by genomic coordinates using samtools
apptainer exec --bind /data "$CONTAINER" \
  samtools sort \
    -@ "$SLURM_CPUS_PER_TASK" \            # Number of CPU threads
    "$BAM_DIR/${SAMPLE}.bam" \             # Input unsorted BAM file
    -o "$SORTED_DIR/${SAMPLE}_sorted.bam"  # Output coordinate-sorted BAM file


# Index the sorted BAM file
# Print status message
echo "Indexing BAM for $SAMPLE"

# Create BAM index file (.bai) for the coordinate-sorted BAM
apptainer exec --bind /data "$CONTAINER" \
  samtools index \
  "$SORTED_DIR/${SAMPLE}_sorted.bam"   # Input sorted BAM file

# Print completion message for this sample
echo "Finished processing $SAMPLE"