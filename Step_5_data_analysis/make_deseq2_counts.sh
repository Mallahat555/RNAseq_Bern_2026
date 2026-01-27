#!/bin/bash
set -euo pipefail  # This line enables strict error handling in Bash.


# set directories as variables
# featureCounts table input directory
FC_DIR=/data/users/apiatkowska/RNAseq/Step_4_count_reads/reads_per_gene_per_sample
# DESeq2-ready tables output
OUT_DIR=/data/users/apiatkowska/RNAseq/Step_5_data_analysis/tables_of_counts_for_DESeq2

# makes output directory if it doesnt exist
mkdir -p "$OUT_DIR"

# =======================
# Clean each sample file
# =======================

# Print a message so the user knows what stage the script is in
echo "Cleaning featureCounts files..."

# Loop over all featureCounts output files
# Each file corresponds to one sample
for FILE in "$FC_DIR"/*_reads_per_gene.txt; do

  # Extract the sample name from the filename
  # Example: FILE = /path/SRR7821949_reads_per_gene.txt
  # SAMPLE = SRR7821949
  SAMPLE=$(basename "$FILE" _reads_per_gene.txt)

  # Remove the first line (featureCounts metadata)
  # Keep only:
  #   column 1 -> Gene ID
  #   column 7 -> read counts for this sample
  #
  # Output is written to a cleaned, per-sample count file
  # Removes the header and keeps Geneid + sample column
  tail -n +2 "$FILE" | cut -f1,7 > "$OUT_DIR/${SAMPLE}_counts.txt"
    
  # ======================
  # Fix the column header:
  # ======================

  # featureCounts uses the full BAM file path as column name
  # This command removes everything up to the last '/'
  sed -i "1s|/data.*/||" "$OUT_DIR/${SAMPLE}_counts.txt"

  # Remove the '_sorted.bam' suffix from the column name
  # Resulting header is just the sample name (e.g. SRR7821949)
  sed -i "1s/_sorted.bam//" "$OUT_DIR/${SAMPLE}_counts.txt"

done

# Print confirmation that all files were processed successfully
echo "Cleaning complete."

# =================
# Extract gene IDs
# =================

# Print a message indicating the current step of the script
echo "Extracting gene IDs..."

# Select the first cleaned per-sample count file
# All *_counts.txt files have identical gene order,
# so we can safely use any one of them
FIRST_SAMPLE=$(ls "$OUT_DIR"/*_counts.txt | head -n 1)

# Extract the first column (Gene IDs) from that file
# and write it to a separate file
# This will be used as the first column of the merged count matrix
cut -f1 "$FIRST_SAMPLE" > "$OUT_DIR/genes.txt"

# ====================
# Extract raw read counts for this sample
# ==================== 

# Print a message indicating the current processing step
echo "Extracting count columns..."

# Initialize an empty Bash array
# This array will store the paths to the per-sample count-only files
COUNT_FILES=()

# Loop over all cleaned per-sample count files
for FILE in "$OUT_DIR"/*_counts.txt; do

  # Extract the sample name from the filename
  # Example:
  # FILE = /path/SRR7821949_counts.txt
  # SAMPLE = SRR7821949
  SAMPLE=$(basename "$FILE" _counts.txt)

  # tail -n +2 "$FILE"   → skip the first line of the featureCounts file (header)
  # cut -f2 "$FILE"      → extract only the second column, which contains the numeric counts
  # > "$OUT_DIR/${SAMPLE}_only.txt" → write the output to a new file in the output directory
  # The resulting file will contain exactly one numeric count per gene,
  # with no header or extra text, ready to merge into the DESeq2 matrix
  tail -n +2 "$FILE" | cut -f2 > "$OUT_DIR/${SAMPLE}_only.txt"

  # Add the path of the count-only file to the COUNT_FILES array
  # This array will later be used to merge all samples into a matrix
  COUNT_FILES+=("$OUT_DIR/${SAMPLE}_only.txt")

done

# ==================
# Merge into matrix
# ==================

# Print a message indicating the current processing step
echo "Merging into count matrix..."

# Merge the gene ID column with all per-sample count columns
# - genes.txt provides the first column (gene identifiers)
# - ${COUNT_FILES[@]} expands to all count-only files, in the same order
# - paste combines columns side-by-side using tabs
#
# The resulting file will contain:
#   Geneid | Sample1 | Sample2 | Sample3 | ...
# but it should NOT include a duplicate header line from genes.txt

# Skip the first line of genes.txt (if it contains a header) using 'tail -n +2'
# Then paste together the gene IDs and all sample count columns
paste <(tail -n +2 "$OUT_DIR/genes.txt") "${COUNT_FILES[@]}" > "$OUT_DIR/counts_no_header.txt"

# Now counts_no_header.txt contains only the data rows (gene IDs + counts)
# The header will be added in the next step

# ==============
# Create header
# ==============

# Print a message indicating the current processing step
echo "Creating header..."

# Initialize the header string with the first column name
# This column will contain gene identifiers
HEADER="Geneid"

# Loop over all cleaned per-sample count files
# The order of this loop defines the column order in the final matrix
for FILE in "$OUT_DIR"/*_counts.txt; do

  # Extract the sample name from the filename
  # Example:
  # FILE = /path/SRR7821949_counts.txt
  # SAMPLE = SRR7821949
  SAMPLE=$(basename "$FILE" _counts.txt)

  # Append the sample name to the header string, separated by a tab
  HEADER="${HEADER}\t${SAMPLE}"

done

# Prepend the header line to the count matrix
# - echo -e prints the header string with interpreted tab characters
# - cat - file concatenates the header (stdin) with the matrix
# - the result is written to the final DESeq2-ready file
echo -e "$HEADER" | \
cat - "$OUT_DIR/counts_no_header.txt" > "$OUT_DIR/counts_all_samples.txt"

# ========
# Cleanup
# ========

# Remove temporary files that are no longer needed:
# - *_only.txt       : the per-sample count-only files
# - genes.txt        : the separate gene ID file
# - counts_no_header.txt : the merged matrix without header
# This keeps the output directory clean and only leaves the final DESeq2 input
rm "$OUT_DIR"/*_only.txt "$OUT_DIR/genes.txt" "$OUT_DIR/counts_no_header.txt"

# ===============
# Final messages
# ===============

# Inform the user that the script has finished successfully
echo "Done!"

# Show the location of the final DESeq2-ready count matrix
echo "Final DESeq2 input:"
echo "$OUT_DIR/counts_all_samples.txt"
