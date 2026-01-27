#!/bin/bash
#SBATCH --job-name=featureCounts_reads_per_gene_per_sample
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --partition=pibu_el8
#SBATCH --array=0-14

#SBATCH --output=/data/users/apiatkowska/RNAseq/Step_4_count_reads/logs/%x_%A_%a.out
#SBATCH --error=/data/users/apiatkowska/RNAseq/Step_4_count_reads/logs/%x_%A_%a.err

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

# Set directories as variables:
# Apptainer container with featureCounts (Subread package)
CONTAINER=/containers/apptainer/subread_2.0.6.sif
# Directory containing coordinate-sorted BAM files (input for featureCounts)
BAM_DIR=/data/users/apiatkowska/RNAseq/Step_3_map_reads/mapped_reads/bam_sorted
# Gene annotation file in GTF format (mouse, GRCm39)
GTF=/data/users/apiatkowska/RNAseq/Step_3_map_reads/reference_genome_and_annotations/Mus_musculus.GRCm39.115.gtf
# Base output directory for read counting results
OUT_DIR=/data/users/apiatkowska/RNAseq/Step_4_count_reads/
# Directory for per-sample gene-level read count tables
SAMPLE_OUT_DIR=$OUT_DIR/reads_per_gene_per_sample
# Directory for standard output and error log files
LOG_DIR=$OUT_DIR/logs

# Create directories if they don't exist
mkdir -p $SAMPLE_OUT_DIR
mkdir -p $LOG_DIR

# Run featureCounts for one sample
# Execute featureCounts inside the Apptainer container
apptainer exec --bind /data "$CONTAINER" featureCounts \
  -T "$SLURM_CPUS_PER_TASK" \        # Number of CPU threads to use
  -p \                               # Specify paired-end reads
  -s 2 \                             # Strand-specific protocol (reverse stranded)
  -a "$GTF" \                        # Gene annotation file in GTF format
  -o "$SAMPLE_OUT_DIR/${SAMPLE}_reads_per_gene.txt" \  # Output gene count table
  "$BAM_DIR/${SAMPLE}_sorted.bam"    # Input coordinate-sorted BAM file

# Confirmation message
echo "Finished counting reads for $SAMPLE"
