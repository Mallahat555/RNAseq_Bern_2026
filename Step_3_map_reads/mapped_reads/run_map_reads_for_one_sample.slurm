#!/bin/bash
#SBATCH --job-name=mapping_hisat2_for_one_sample
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --time=04:00:00
#SBATCH --partition=pibu_el8

#SBATCH --output=/data/users/apiatkowska/RNAseq/Step_3_map_reads/mapped_reads/logs/%x_%j.out
#SBATCH --error=/data/users/apiatkowska/RNAseq/Step_3_map_reads/mapped_reads/logs/%x_%j.err

#SBATCH --mail-user=agnieszka.piatkowska@students.unibe.ch
#SBATCH --mail-type=BEGIN,END

# Sample name
SAMPLE=SRR7821949

# Directories
CONTAINER=/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif
INDEX_DIR=/data/users/apiatkowska/RNAseq/Step_3_map_reads/hisat2_indexing/index/mm39
FASTQ_DIR=/data/users/apiatkowska/RNAseq/Step_2.2_trimmed_fastQ/trimmed_fastQ
OUT_DIR=/data/users/apiatkowska/RNAseq/Step_3_map_reads/mapped_reads

# Create directory structure
mkdir -p $OUT_DIR/sam
mkdir -p $OUT_DIR/bam
mkdir -p $OUT_DIR/bam_sorted
mkdir -p $OUT_DIR/logs

# Define FASTQ files for the sample
FASTQ_R1=${FASTQ_DIR}/${SAMPLE}_1_trimmed.fastq.gz
FASTQ_R2=${FASTQ_DIR}/${SAMPLE}_2_trimmed.fastq.gz

# Map reads with HISAT2 (strand-specific RF)
apptainer exec --bind /data/ $CONTAINER hisat2 \
  -x $INDEX_DIR \
  -1 $FASTQ_R1 \
  -2 $FASTQ_R2 \
  --rna-strandness RF \
  -p $SLURM_CPUS_PER_TASK \
  -S $OUT_DIR/sam/${SAMPLE}.sam

# Convert SAM → BAM
apptainer exec --bind /data/ $CONTAINER samtools view -bS $OUT_DIR/sam/${SAMPLE}.sam > $OUT_DIR/bam/${SAMPLE}.bam

# Sort BAM by genomic coordinates
apptainer exec --bind /data/ $CONTAINER samtools sort $OUT_DIR/bam/${SAMPLE}.bam -o $OUT_DIR/bam_sorted/${SAMPLE}_sorted.bam

# Index the sorted BAM
apptainer exec --bind /data/ $CONTAINER samtools index $OUT_DIR/bam_sorted/${SAMPLE}_sorted.bam


