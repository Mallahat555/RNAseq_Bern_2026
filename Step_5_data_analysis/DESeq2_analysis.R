# Load DESeq2 and other useful packages
library(DESeq2)
library(ggplot2)
library(pheatmap)

# verify if libraries are loaded
sessionInfo()

# Define the path to your counts table and set it into the variable
counts_file <- "C:/Users/agnie/Desktop/MSc Bioinformatics/MSc Bioinformatics Courses/RNA-seq/Project_RNA_seq/Step_5_data_analysis/tables_of_counts_for_DESeq2/counts_all_samples.txt"

# Read the counts table from counts_file:
# read.delim() Each line in your file becomes a row in a data frame.
# row.names = 1 use the first column as row names instead of creating a separate numeric index
# check.names = FALSE keeps the sample names exactly as they appear in file
counts <- read.delim(counts_file, row.names=1, check.names=FALSE)

# Check the first few rows - just a control if all is ok
head(counts)

#----------------------------------------------------------------------------
# Create the sample metadata and set up the DESeq2 dataset.
#----------------------------------------------------------------------------

# Sample names (columns in counts)
samples <- colnames(counts)

# Define the group for each sample in the same order as columns
group <- c(
  "Blood_WT_Case", "Blood_WT_Case", "Blood_WT_Case", "Blood_WT_Case", "Blood_WT_Case",
  "Blood_DKO_Case", "Blood_DKO_Case", "Blood_DKO_Case", "Blood_DKO_Case",
  "Blood_WT_Control", "Blood_WT_Control", "Blood_WT_Control",
  "Blood_DKO_Control", "Blood_DKO_Control", "Blood_DKO_Control"
)

# Create metadata dataframe
colData <- data.frame(row.names = samples,
                      group = factor(group))

# Check the metadata
colData

#-------------------------------
# Create DESeqDataSet object
#-------------------------------

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = counts,      # raw counts table
  colData = colData,       # sample metadata
  design = ~ group         # experimental design formula
)

# Check the DESeqDataSet object
dds

#---------------------------
#Pre-filter low-count genes
#---------------------------

# Keep genes with at least 10 reads total across all samples
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Check dimensions after filtering
dim(dds)

#--------------------------------------------------------------
# Histograms Gene-wise total counts before and after filtering
#--------------------------------------------------------------

# Counts BEFORE filtering
counts_before <- counts
total_counts_before <- rowSums(counts_before)

# Counts AFTER filtering
counts_after <- counts(dds)
total_counts_after <- rowSums(counts_after)

# Compute log10 counts
log_counts_before <- log10(total_counts_before + 1)
log_counts_after  <- log10(total_counts_after + 1)

# Decide the break in y-axis
break_low  <- 0       # lower part
break_high <- 2500    # upper part of first plot
break_top  <- 40000   # top bar

# First plot: normal range
hist(log_counts_before, breaks=50,
     col=rgb(0,0,1,0.4), border="white",
     xlab="log10 of total reads per gene", ylab="frequency distribution",
     main="Effect of filtering on reads counts",
     ylim=c(break_low, break_high))
hist(log_counts_after, breaks=50,
     col=rgb(0,1,0,0.4), border="white",
     add=TRUE)

# Add legend
legend("topright",
       legend=c("before filtering", "after filtering"),
       fill=c(rgb(0,0,1,0.4), rgb(0,1,0,0.4)),
       border="white")

#-----------------------------------------
# Run DESeq normalization & model fitting
#-----------------------------------------

dds_raw <- dds        # keep original filtered counts (in case I screw something)
dds <- DESeq(dds_raw) # run DESeq on a copy

#------------------------------------------
# Inspect counts and transformations
#------------------------------------------

# Variance stabilizing transformation
vsd <- vst(dds, blind=TRUE)

# Or rlog transformation
rld <- rlog(dds, blind=TRUE)

#----------------
# Quality Checks
#----------------

#------------------------------------------
# PCA plot 
#------------------------------------------

# Compute PCA and extract data (returns coordinates and variance)
pca_data <- plotPCA(vsd, intgroup="group", returnData = TRUE)

# Extract percentage variance explained by PC1 and PC2
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Remove "Blood_" prefix for legend labels only
legend_labels <- gsub("^Blood_", "", levels(pca_data$group))

# Build PCA plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4) +
  
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  
  # Adjusted "Blood" label in top-left
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = "Blood",
    hjust = -0.12,  # horizontal position of the label
    vjust = 1.5,    # vertical position of the label 
    size = 6,
    fontface = "bold"
  ) +
  
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12)
  ) +
  
  coord_fixed(ratio = 1) +
  
  # Use modified labels in legend only
  scale_color_discrete(labels = legend_labels)

# Display PCA plot
pca_plot

#--------------------------
# Sample distance Heat map
#--------------------------
colData(vsd)

library(RColorBrewer)
library(pheatmap)

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

sampleLabels <- paste(
  rownames(colData(vsd)),
  colData(vsd)$group,
  sep = " | "
)

rownames(sampleDistMatrix) <- sampleLabels
colnames(sampleDistMatrix) <- sampleLabels

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  col = colors,
  main = "Sample-to-sample distance heatmap"
)

#-----------------------------
# Checking groups and contrast
#-----------------------------

# check groups
design = ~ group
levels(dds$group)

# what pairwise comparisons are avaiable
resultsNames(dds)

#----------
# Contrasts
#----------

# 1. WT Case vs WT Control (up in WT Case)
res_WT_case_vs_WT_control <- results(
  dds,
  contrast = c("group", "Blood_WT_Case", "Blood_WT_Control")
)

# 2. DKO Case vs DKO Control (up in DKO Case)
res_DKO_case_vs_DKO_control <- results(
  dds,
  contrast = c("group", "Blood_DKO_Case", "Blood_DKO_Control")
)

# 3. DKO Case vs WT Case (up in DKO Case)
res_DKO_case_vs_WT_case <- results(
  dds,
  contrast = c("group", "Blood_DKO_Case", "Blood_WT_Case")
)


#--------------------
# Look at the results
#--------------------

#---------------------------------
# 3. DKO Case vs WT Case (up in DKO Case)
#---------------------------------

# View first few genes
head(res_DKO_case_vs_WT_case)

# Open in RStudio table viewer
View(res_DKO_case_vs_WT_case)

# Column names
colnames(res_DKO_case_vs_WT_case)

#-------------------------------
# Filter significant genes
# padj < 0.05 and |log2FC| > 1
#-------------------------------
res_sig <- res_DKO_case_vs_WT_case[
  !is.na(res_DKO_case_vs_WT_case$padj) &
    res_DKO_case_vs_WT_case$padj < 0.05 &
    abs(res_DKO_case_vs_WT_case$log2FoldChange) > 1,
]

# Number of significant genes
nrow(res_sig)

# Up- vs down-regulated genes
up_genes   <- sum(res_sig$log2FoldChange > 0)  # higher in DKO_case
down_genes <- sum(res_sig$log2FoldChange < 0)  # higher in WT_case

up_genes
down_genes

#-------------------------------
# Summary (all genes, for context)
#-------------------------------
summary(res_DKO_case_vs_WT_case, alpha = 0.05)

#-------------------------------
# MA plot
#-------------------------------
plotMA(
  res_DKO_case_vs_WT_case,
  xlab = "Mean normalized expression",
  ylab = "Log2 fold change",
  main = "DKO Case vs WT Case"
)

#-------------------------------
# Save significant genes to CSV
#-------------------------------
write.csv(
  as.data.frame(res_sig),
  file = "DKO_case_vs_WT_case_DESeq2_sig_genes.csv"
)
#-------------------------------------
# 2. DKO Case vs DKO Control (up in DKO Case)
#-------------------------------------

# View first few genes
head(res_DKO_case_vs_DKO_control)

# Open in RStudio table viewer
View(res_DKO_case_vs_DKO_control)

# Column names
colnames(res_DKO_case_vs_DKO_control)

#-------------------------------
# Filter significant genes
# padj < 0.05 and |log2FC| > 1
#-------------------------------
res_sig_DKO_case_vs_DKO_control <- res_DKO_case_vs_DKO_control[
  !is.na(res_DKO_case_vs_DKO_control$padj) &
    res_DKO_case_vs_DKO_control$padj < 0.05 &
    abs(res_DKO_case_vs_DKO_control$log2FoldChange) > 1,
]

# Number of significant genes
nrow(res_sig_DKO_case_vs_DKO_control)

# Up- vs down-regulated genes
up_genes_DKO   <- sum(res_sig_DKO_case_vs_DKO_control$log2FoldChange > 0)  # higher in DKO_Case
down_genes_DKO <- sum(res_sig_DKO_case_vs_DKO_control$log2FoldChange < 0)  # higher in DKO_Control

up_genes_DKO
down_genes_DKO

#-------------------------------
# Summary (all genes, for context)
#-------------------------------
summary(res_DKO_case_vs_DKO_control, alpha = 0.05)

#-------------------------------
# MA plot
#-------------------------------
plotMA(
  res_DKO_case_vs_DKO_control,
  xlab = "Mean normalized expression",
  ylab = "Log2 fold change",
  main = "DKO Case vs DKO Control"
)

#-------------------------------
# Save significant genes to CSV
#-------------------------------
write.csv(
  as.data.frame(res_sig_DKO_case_vs_DKO_control),
  file = "DKO_case_vs_DKO_control_DESeq2_sig_genes.csv"
)

#-------------------------------------
# 1. WT Case vs WT Control (up in WT Case)
#-------------------------------------

# View first few genes
head(res_WT_case_vs_WT_control)

# Open in RStudio table viewer
View(res_WT_case_vs_WT_control)

# Column names
colnames(res_WT_case_vs_WT_control)

#-------------------------------
# Filter significant genes
# padj < 0.05 and |log2FC| > 1
#-------------------------------
res_sig_WT_case_vs_WT_control <- res_WT_case_vs_WT_control[
  !is.na(res_WT_case_vs_WT_control$padj) &
    res_WT_case_vs_WT_control$padj < 0.05 &
    abs(res_WT_case_vs_WT_control$log2FoldChange) > 1,
]

# Number of significant genes
nrow(res_sig_WT_case_vs_WT_control)

# Up- vs down-regulated genes
up_genes_WT   <- sum(res_sig_WT_case_vs_WT_control$log2FoldChange > 0)  # higher in WT_Case
down_genes_WT <- sum(res_sig_WT_case_vs_WT_control$log2FoldChange < 0)  # higher in WT_Control

up_genes_WT
down_genes_WT

#-------------------------------
# Summary (all genes, for context)
#-------------------------------
summary(res_WT_case_vs_WT_control, alpha = 0.05)

#-------------------------------
# MA plot
#-------------------------------
plotMA(
  res_WT_case_vs_WT_control,
  xlab = "Mean normalized expression",
  ylab = "Log2 fold change",
  main = "WT Case vs WT Control"
)

#-------------------------------
# Save significant genes to CSV
#-------------------------------
write.csv(
  as.data.frame(res_sig_WT_case_vs_WT_control),
  file = "WT_case_vs_WT_control_DESeq2_sig_genes.csv"
)
#----------------------------
# Functions for Venn diagrams
#----------------------------
install.packages("VennDiagram")
install.packages("gridExtra")  # usually already installed, but safe
library(VennDiagram)

#-------------------------------
# Significant DE genes (|log2FC|>1 & padj<0.05)
#-------------------------------

genes_WT_case_vs_WT_control <- rownames(res_WT_case_vs_WT_control[
  !is.na(res_WT_case_vs_WT_control$padj) &
    res_WT_case_vs_WT_control$padj < 0.05 &
    abs(res_WT_case_vs_WT_control$log2FoldChange) > 1,
])

genes_DKO_case_vs_DKO_control <- rownames(res_DKO_case_vs_DKO_control[
  !is.na(res_DKO_case_vs_DKO_control$padj) &
    res_DKO_case_vs_DKO_control$padj < 0.05 &
    abs(res_DKO_case_vs_DKO_control$log2FoldChange) > 1,
])

genes_DKO_case_vs_WT_case <- rownames(res_DKO_case_vs_WT_case[
  !is.na(res_DKO_case_vs_WT_case$padj) &
    res_DKO_case_vs_WT_case$padj < 0.05 &
    abs(res_DKO_case_vs_WT_case$log2FoldChange) > 1,
])

#-------------------------------
# Venn diagram for all significant DE genes
#-------------------------------
venn.plot <- venn.diagram(
  x = list(
    "WT Case vs WT Control"   = genes_WT_case_vs_WT_control,
    "DKO Case vs DKO Control" = genes_DKO_case_vs_DKO_control,
    "DKO Case vs WT Case"     = genes_DKO_case_vs_WT_case
  ),
  filename = NULL,
  fill = c("skyblue", "pink", "lightgreen"),
  alpha = 0.6,
  cex = 1.2,
  cat.cex = 1.2
)

grid::grid.draw(venn.plot)

#-------------------------------
# UP-regulated genes (log2FC > 1)
#-------------------------------
up_WT_case_vs_WT_control <- rownames(res_WT_case_vs_WT_control[
  !is.na(res_WT_case_vs_WT_control$padj) &
    res_WT_case_vs_WT_control$padj < 0.05 &
    res_WT_case_vs_WT_control$log2FoldChange > 1,
])

up_DKO_case_vs_DKO_control <- rownames(res_DKO_case_vs_DKO_control[
  !is.na(res_DKO_case_vs_DKO_control$padj) &
    res_DKO_case_vs_DKO_control$padj < 0.05 &
    res_DKO_case_vs_DKO_control$log2FoldChange > 1,
])

up_DKO_case_vs_WT_case <- rownames(res_DKO_case_vs_WT_case[
  !is.na(res_DKO_case_vs_WT_case$padj) &
    res_DKO_case_vs_WT_case$padj < 0.05 &
    res_DKO_case_vs_WT_case$log2FoldChange > 1,
])

venn_up <- venn.diagram(
  x = list(
    "WT Case vs WT Control"   = up_WT_case_vs_WT_control,
    "DKO Case vs DKO Control" = up_DKO_case_vs_DKO_control,
    "DKO Case vs WT Case"     = up_DKO_case_vs_WT_case
  ),
  filename = NULL,
  fill = c("skyblue", "pink", "lightgreen"),
  alpha = 0.6,
  cex = 3,
  cat.cex = 1.8,
  main = "Up-regulated genes",
  main.cex = 2,
  cat.pos = c(-20, 20, 180),
  cat.dist = c(0.05, 0.05, 0.05)
)

grid::grid.newpage()
grid::grid.draw(venn_up)

#-------------------------------
# DOWN-regulated genes (log2FC < -1)
#-------------------------------
down_WT_case_vs_WT_control <- rownames(res_WT_case_vs_WT_control[
  !is.na(res_WT_case_vs_WT_control$padj) &
    res_WT_case_vs_WT_control$padj < 0.05 &
    res_WT_case_vs_WT_control$log2FoldChange < -1,
])

down_DKO_case_vs_DKO_control <- rownames(res_DKO_case_vs_DKO_control[
  !is.na(res_DKO_case_vs_DKO_control$padj) &
    res_DKO_case_vs_DKO_control$padj < 0.05 &
    res_DKO_case_vs_DKO_control$log2FoldChange < -1,
])

down_DKO_case_vs_WT_case <- rownames(res_DKO_case_vs_WT_case[
  !is.na(res_DKO_case_vs_WT_case$padj) &
    res_DKO_case_vs_WT_case$padj < 0.05 &
    res_DKO_case_vs_WT_case$log2FoldChange < -1,
])

venn_down <- venn.diagram(
  x = list(
    "WT Case vs WT Control"   = down_WT_case_vs_WT_control,
    "DKO Case vs DKO Control" = down_DKO_case_vs_DKO_control,
    "DKO Case vs WT Case"     = down_DKO_case_vs_WT_case
  ),
  filename = NULL,
  fill = c("skyblue", "pink", "lightgreen"),
  alpha = 0.6,
  cex = 3,
  cat.cex = 1.8,
  main = "Down-regulated genes",
  main.cex = 2,
  cat.pos = c(-20, 20, 180),
  cat.dist = c(0.05, 0.05, 0.05)
)

grid::grid.newpage()
grid::grid.draw(venn_down)

#-------------------------------------------------------------
#Select 2-3 genes of particular interest from the publication
#-------------------------------------------------------------

# Install the Bioconductor package for mouse gene annotation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Mm.eg.db")

# Make sure AnnotationDbi is installed
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  BiocManager::install("AnnotationDbi")
}

# Load the package
library(org.Mm.eg.db)

# Load AnnotationDbi
library(AnnotationDbi)

#-------------------------------
# Selected genes of interest
#-------------------------------
genes_of_interest <- c("Tap1", "Irf7", "Oas1a")

gene_map <- mapIds(
  org.Mm.eg.db,
  keys = genes_of_interest,
  column = "ENSEMBL",
  keytype = "SYMBOL",
  multiVals = "first"
)

# Extract Ensembl IDs as a plain character vector
genes_of_interest_ensembl <- as.character(gene_map)

genes_of_interest_ensembl

#-------------------------------
# Check presence in DESeq2 results
#-------------------------------
genes_present <- list(
  WT_Case_vs_WT_Control   = genes_of_interest_ensembl %in% rownames(res_WT_case_vs_WT_control),
  DKO_Case_vs_DKO_Control = genes_of_interest_ensembl %in% rownames(res_DKO_case_vs_DKO_control),
  DKO_Case_vs_WT_Case     = genes_of_interest_ensembl %in% rownames(res_DKO_case_vs_WT_case)
)

# Name the TRUE/FALSE vector by gene symbols for clarity
genes_present <- lapply(genes_present, function(x) { names(x) <- genes_of_interest; x })

genes_present
#-----------------------------------------
# Gene labels: Ensembl ID -> Gene symbol
#-----------------------------------------
label_map <- data.frame(
  ensembl_id = c(
    "ENSMUSG00000037321",  # Tap1
    "ENSMUSG00000025498",  # Irf7
    "ENSMUSG00000052776"   # Oas1
  ),
  symbol = c("Tap1", "Irf7", "Oas1"),
  stringsAsFactors = FALSE
)

genes_of_interest_ensembl <- label_map$ensembl_id

#---------------------------------
# Function to make a volcano plot
#---------------------------------
install.packages("ggrepel", repos="https://cloud.r-project.org")
library(ggrepel)  # for non-overlapping labels

plot_volcano_ensembl <- function(res,
                                 title,
                                 genes_of_interest,
                                 label_df,
                                 padj_cutoff = 0.05,
                                 lfc_cutoff = 1) {
  
  df <- as.data.frame(res)
  df$ensembl_id <- rownames(df)
  
  # Replace NA padj for plotting
  df$padj_plot <- ifelse(is.na(df$padj), 1, df$padj)
  
  # Base gene status
  df$status <- "Not significant"
  df$status[df$padj < padj_cutoff & df$log2FoldChange >  lfc_cutoff]  <- "Up-regulated"
  df$status[df$padj < padj_cutoff & df$log2FoldChange < -lfc_cutoff] <- "Down-regulated"
  
  # Override for genes of interest
  df$status[df$ensembl_id %in% genes_of_interest] <- "Gene of interest"
  
  df$status <- factor(
    df$status,
    levels = c("Up-regulated",
               "Down-regulated",
               "Not significant",
               "Gene of interest")
  )
  
  # Split for layering
  df_bg <- df[df$status != "Gene of interest", ]
  df_fg <- df[df$status == "Gene of interest", ]
  
  # Merge symbols for labels
  df_labels <- merge(
    df_fg,
    label_df,
    by = "ensembl_id",
    all.x = TRUE
  )
  
  ggplot() +
    ## Background genes
    geom_point(
      data = df_bg,
      aes(x = log2FoldChange,
          y = -log10(padj_plot),
          color = status),
      size = 2,
      alpha = 0.6
    ) +
    
    ## Genes of interest (yellow, on top)
    geom_point(
      data = df_fg,
      aes(x = log2FoldChange,
          y = -log10(padj_plot),
          color = status),
      size = 2
    ) +
    
    ## Gene symbol labels only
    geom_text_repel(
      data = df_labels,
      aes(x = log2FoldChange,
          y = -log10(padj_plot),
          label = symbol),
      size = 4,
      box.padding = 0.5,
      max.overlaps = Inf
    ) +
    
    ## Threshold lines
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
               linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff),
               linetype = "dashed") +
    
    ## Colors + legend
    scale_color_manual(
      values = c(
        "Up-regulated"     = "firebrick",
        "Down-regulated"   = "royalblue",
        "Not significant"  = "grey70",
        "Gene of interest" = "gold"
      ),
      drop = FALSE
    ) +
    
    labs(
      title = title,
      x = "Log2 fold change",
      y = "-log10(adjusted p-value)",
      color = "Gene category"
    ) +
    
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      # Legend at top-left
      legend.position = c(0.05, 0.95),       # top-left corner
      legend.justification = c("left", "top"),
      legend.background = element_rect(fill = alpha("white", 0.6)) # optional semi-transparent background
    )
}

#-----------------------------------------
# Create volcano plots for each contrast
#-----------------------------------------

# 1. WT Case vs WT Control
volcano_WT <- plot_volcano_ensembl(
  res = res_WT_case_vs_WT_control,
  title = "WT case vs WT control",
  genes_of_interest = genes_of_interest_ensembl,
  label_df = label_map
)
volcano_WT

# 2. DKO Case vs DKO Control
volcano_DKO <- plot_volcano_ensembl(
  res = res_DKO_case_vs_DKO_control,
  title = "DKO case vs DKO control",
  genes_of_interest = genes_of_interest_ensembl,
  label_df = label_map
)
volcano_DKO

# 3. DKO Case vs WT Case
volcano_DKOvsWT <- plot_volcano_ensembl(
  res = res_DKO_case_vs_WT_case,
  title = "DKO case vs WT case",
  genes_of_interest = genes_of_interest_ensembl,
  label_df = label_map
)
volcano_DKOvsWT


#-------------------------------------------------
# Extract normalised counts for genes of interest
#-------------------------------------------------

# Normalized counts for all genes
norm_counts <- counts(dds, normalized = TRUE)

# Subset for your genes of interest
norm_counts_goi <- norm_counts[gene_map, ]

# Inspect
norm_counts_goi

#----------------------
# Convert for plotting
#----------------------

# Make sure tidyr is loaded
library(tidyr)
library(dplyr)  # for convenient data manipulation

# Convert normalized counts to a tidy long format
norm_counts_long <- as.data.frame(norm_counts_goi) %>%
  tibble::rownames_to_column(var = "ensembl_id") %>%  # move rownames to a column
  pivot_longer(
    cols = -ensembl_id,
    names_to = "sample",
    values_to = "normalized_counts"
  )

# Add gene symbol using your label map
norm_counts_long$gene <- label_map$symbol[match(norm_counts_long$ensembl_id, label_map$ensembl_id)]

# Add group info from colData
norm_counts_long$group <- colData[norm_counts_long$sample, "group"]

# Check
head(norm_counts_long)

#------------------------------------------------------------
# Plot normalised expressions for selected genes of interest
#------------------------------------------------------------

# Create a clean group name column for plotting
norm_counts_long <- norm_counts_long %>%
  mutate(group_clean = gsub("^Blood_", "", group),   # remove Blood_ prefix
         group_clean = gsub("_", " ", group_clean))  # replace underscores with space

# Plot normalized counts
ggplot(norm_counts_long, aes(x = group_clean, y = normalized_counts, color = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(width = 0.2, size = 2) +
  facet_wrap(~ gene, scales = "free_y") +
  theme_bw() +
  labs(
    y = "Normalized counts",
    x = "Group",
    title = "Expression of selected genes"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

#-------------------------------------------------
# Extract normalized counts for genes of interest
#-------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
norm_counts <- counts(dds, normalized = TRUE)
norm_counts_goi <- norm_counts[gene_map, ]

norm_counts_long <- as.data.frame(norm_counts_goi) %>%
  tibble::rownames_to_column(var = "ensembl_id") %>%
  pivot_longer(
    cols = -ensembl_id,
    names_to = "sample",
    values_to = "normalized_counts"
  )

norm_counts_long$gene <- label_map$symbol[match(norm_counts_long$ensembl_id, label_map$ensembl_id)]
norm_counts_long$group <- colData[norm_counts_long$sample, "group"]
norm_counts_long <- norm_counts_long %>%
  mutate(group_clean = gsub("^Blood_", "", group),
         group_clean = gsub("_", " ", group_clean),
         group_clean = factor(group_clean, levels = c(
           "WT Control", "WT Case", "DKO Control", "DKO Case"
         )))

#--------------------------------------------
# Define comparisons
#--------------------------------------------
comparisons <- list(
  c("WT Control", "WT Case"),
  c("DKO Control", "DKO Case"),
  c("WT Case", "DKO Case")
)

#--------------------------------------------
# Compute Wilcoxon p-values and positions per gene
#--------------------------------------------
stat_results <- do.call(rbind, lapply(unique(norm_counts_long$gene), function(g) {
  gene_data <- norm_counts_long %>% filter(gene == g)
  
  # max y for this gene
  y_max <- max(gene_data$normalized_counts)
  
  # offsets for multiple comparisons per gene
  offsets <- seq(0.05, 0.15, length.out = length(comparisons))
  
  # create data frame for all comparisons of this gene
  df <- data.frame(
    gene = g,
    comp1 = sapply(comparisons, "[", 1),
    comp2 = sapply(comparisons, "[", 2),
    y_bracket = y_max * (1 + offsets),         # bracket line
    y_signif = y_max * (1 + offsets + 0.02),  # asterisk/ns
    y_pval   = y_max * (1 + offsets + 0.04)   # p-value above asterisk
  )
  
  # compute p-values and significance
  df <- df %>%
    rowwise() %>%
    mutate(
      values1 = list(gene_data$normalized_counts[gene_data$group_clean == comp1]),
      values2 = list(gene_data$normalized_counts[gene_data$group_clean == comp2]),
      p_value = wilcox.test(unlist(values1), unlist(values2))$p.value,
      signif = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    ) %>%
    ungroup()
  
  return(df)
}))

#--------------------------------------------
# Plot
#--------------------------------------------
ggplot(norm_counts_long, aes(x = group_clean, y = normalized_counts, color = group_clean)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(width = 0.2, size = 2) +
  facet_wrap(~ gene, scales = "free_y") +
  # Bracket lines
  geom_segment(data = stat_results,
               aes(x = comp1, xend = comp2, y = y_bracket, yend = y_bracket),
               inherit.aes = FALSE, color = "black") +
  geom_segment(data = stat_results,
               aes(x = comp1, xend = comp1, y = y_bracket, yend = y_bracket - 0.01 * y_bracket),
               inherit.aes = FALSE, color = "black") +
  geom_segment(data = stat_results,
               aes(x = comp2, xend = comp2, y = y_bracket, yend = y_bracket - 0.01 * y_bracket),
               inherit.aes = FALSE, color = "black") +
  # Asterisk/ns
  geom_text(data = stat_results,
            aes(x = (as.numeric(factor(comp1, levels = levels(norm_counts_long$group_clean))) +
                       as.numeric(factor(comp2, levels = levels(norm_counts_long$group_clean))))/2,
                y = y_signif,
                label = signif),
            inherit.aes = FALSE, size = 3, color = "black") +
  # p-values
  geom_text(data = stat_results,
            aes(x = (as.numeric(factor(comp1, levels = levels(norm_counts_long$group_clean))) +
                       as.numeric(factor(comp2, levels = levels(norm_counts_long$group_clean))))/2,
                y = y_pval,
                label = paste0("p = ", signif(p_value,3))),
            inherit.aes = FALSE, size = 3, color = "black") +
  theme_bw() +
  labs(
    y = "Normalized counts",
    x = "Group",
    title = "Expression of selected genes with Wilcoxon test significance"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

    
#==================================================
# Overrepresentation Analysis with clusterProfiler
#==================================================

#----------------------------
# 1. Install & load packages
#----------------------------
if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Mm.eg.db", quietly = TRUE))
  BiocManager::install("org.Mm.eg.db")

library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(stringr)
library(dplyr)

#====================================================
# Overrepresentation Analysis: WT case vs WT control
#====================================================
#----------------------------
# 1. Prepare gene lists
#----------------------------

# DE genes for WT case vs WT control (flip contrast)
DE_genes_WT <- rownames(res_sig_WT_case_vs_WT_control)
all_genes_WT <- rownames(res_WT_case_vs_WT_control)

# Remove any NA IDs
DE_genes_WT <- DE_genes_WT[!is.na(DE_genes_WT)]
all_genes_WT <- all_genes_WT[!is.na(all_genes_WT)]

cat("Number of DE genes (WT case vs control):", length(DE_genes_WT), "\n")
cat("Number of all genes (WT case vs control):", length(all_genes_WT), "\n")
head(DE_genes_WT)
head(all_genes_WT)

#------------------------------------
# 2. Run overrepresentation analysis
#------------------------------------
ego_WT <- enrichGO(
  gene          = DE_genes_WT,
  universe      = all_genes_WT,
  OrgDb         = org.Mm.eg.db,
  ont           = "ALL",
  keyType       = "ENSEMBL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

#----------------------------------
# 3. Convert results to data frame
#----------------------------------
ego_WT_df <- as.data.frame(ego_WT)

# Calculate numeric GeneRatio and BgRatio
ego_WT_df$GeneRatioNum <- sapply(ego_WT_df$GeneRatio, function(x) {
  parts <- unlist(strsplit(x, "/"))
  as.numeric(parts[1]) / as.numeric(parts[2])
})
ego_WT_df$BgRatioNum <- sapply(ego_WT_df$BgRatio, function(x) {
  parts <- unlist(strsplit(x, "/"))
  as.numeric(parts[1]) / as.numeric(parts[2])
})

#------------------------------------
# 4. Sort and select top 15 GO terms
#------------------------------------
ego_top_WT <- ego_WT@result %>%
  arrange(desc(Count)) %>%    # largest gene counts first
  slice(1:15)

# Wrap long GO term descriptions
ego_top_WT$Description_wrapped <- str_wrap(ego_top_WT$Description, width = 40)

# Set factor levels so largest dots appear at top
ego_top_WT$Description_wrapped <- factor(
  ego_top_WT$Description_wrapped,
  levels = rev(ego_top_WT$Description_wrapped)
)

#------------------------------------
# 5. Visualize top 15 GO terms
#------------------------------------
ggplot(ego_top_WT, aes(x = GeneRatio, y = Description_wrapped)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_continuous(low = "red", high = "blue", trans = "reverse") +
  scale_size(range = c(3, 8)) +
  theme_bw() +
  labs(
    title = "Top 15 enriched GO terms (WT case vs WT control)",
    x = "Gene Ratio",
    y = NULL,
    color = "Adjusted p-value",
    size = "Gene count"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10, hjust = 1),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # vertical x-axis labels
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )

#=======================================================
# Overrepresentation Analysis: DKO case vs DKO control
#=======================================================

#----------------------------
# 1. Prepare gene lists
#----------------------------
# DE genes for DKO case vs DKO control comparison
DE_genes_DKO <- rownames(res_sig_DKO_case_vs_DKO_control)
all_genes_DKO <- rownames(res_DKO_case_vs_DKO_control)

# Remove any NA IDs
DE_genes_DKO <- DE_genes_DKO[!is.na(DE_genes_DKO)]
all_genes_DKO <- all_genes_DKO[!is.na(all_genes_DKO)]

# Optional checks
cat("Number of DE genes (DKO case vs control):", length(DE_genes_DKO), "\n")
cat("Number of all genes (DKO case vs control):", length(all_genes_DKO), "\n")
head(DE_genes_DKO)
head(all_genes_DKO)

#------------------------------------
# 2. Run overrepresentation analysis
#------------------------------------
ego_DKO <- enrichGO(
  gene          = DE_genes_DKO,         # DE genes
  universe      = all_genes_DKO,        # background genes
  OrgDb         = org.Mm.eg.db,         # mouse annotation database
  ont           = "ALL",                 # GO category
  keyType       = "ENSEMBL",            # our IDs are Ensembl
  pAdjustMethod = "BH",                 # multiple testing correction
  pvalueCutoff  = 0.05,                 # significance threshold
  qvalueCutoff  = 0.05,                 # FDR threshold
  readable      = TRUE                  # converts Ensembl IDs to gene symbols
)

#----------------------------------
# 3. Convert results to data frame
#----------------------------------
ego_DKO_df <- as.data.frame(ego_DKO)

# Calculate numeric GeneRatio and BgRatio
ego_DKO_df$GeneRatioNum <- sapply(ego_DKO_df$GeneRatio, function(x) {
  parts <- unlist(strsplit(x, "/"))
  as.numeric(parts[1]) / as.numeric(parts[2])
})
ego_DKO_df$BgRatioNum <- sapply(ego_DKO_df$BgRatio, function(x) {
  parts <- unlist(strsplit(x, "/"))
  as.numeric(parts[1]) / as.numeric(parts[2])
})

#------------------------------------
# 4. Sort and report top 15 GO terms
#------------------------------------
top10_GO_DKO <- ego_DKO_df[order(ego_DKO_df$p.adjust), ][1:15,
                                                         c("ID", "Description", "GeneRatio", "GeneRatioNum",
                                                           "BgRatio", "BgRatioNum", "p.adjust")]
top10_GO_DKO

#----------------------------
# 5. Visualize results
#----------------------------
library(ggplot2)
library(stringr)
library(dplyr)

# Prepare top 15 terms
ego_top <- ego_DKO@result %>%
  arrange(desc(Count)) %>%      # largest gene counts first
  slice(1:15)

# Wrap long GO term descriptions
ego_top$Description_wrapped <- str_wrap(ego_top$Description, width = 40)

# Set factor levels so largest dots appear at top
ego_top$Description_wrapped <- factor(
  ego_top$Description_wrapped,
  levels = rev(ego_top$Description_wrapped)
)

# Plot
ggplot(ego_top, aes(x = GeneRatio, y = Description_wrapped)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_continuous(low = "red", high = "blue", trans = "reverse") +
  scale_size(range = c(3, 8)) +
  theme_bw() +
  labs(
    title = "Top 15 enriched GO terms (DKO case vs DKO control)",
    x = "Gene Ratio",
    y = NULL,
    color = "Adjusted p-value",
    size = "Gene count"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10, hjust = 1),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # vertical x-axis labels
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )

#====================================================
# Overrepresentation Analysis: DKO case vs WT case
#====================================================

#----------------------------
# 1. Prepare gene lists
#----------------------------
DE_genes_DKOvsWT <- rownames(res_sig)            # significant DE genes
all_genes_DKOvsWT <- rownames(res_DKO_case_vs_WT_case)  # all genes in contrast

# Remove any NA IDs
DE_genes_DKOvsWT <- DE_genes_DKOvsWT[!is.na(DE_genes_DKOvsWT)]
all_genes_DKOvsWT <- all_genes_DKOvsWT[!is.na(all_genes_DKOvsWT)]

cat("Number of DE genes (DKO case vs WT case):", length(DE_genes_DKOvsWT), "\n")
cat("Number of all genes (DKO case vs WT case):", length(all_genes_DKOvsWT), "\n")
head(DE_genes_DKOvsWT)
head(all_genes_DKOvsWT)

#------------------------------------
# 2. Run overrepresentation analysis
#------------------------------------
ego_DKOvsWT <- enrichGO(
  gene          = DE_genes_DKOvsWT,
  universe      = all_genes_DKOvsWT,
  OrgDb         = org.Mm.eg.db,
  ont           = "ALL",
  keyType       = "ENSEMBL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

#----------------------------------
# 3. Convert results to data frame
#----------------------------------
ego_DKOvsWT_df <- as.data.frame(ego_DKOvsWT)

# Calculate numeric GeneRatio and BgRatio
ego_DKOvsWT_df$GeneRatioNum <- sapply(ego_DKOvsWT_df$GeneRatio, function(x) {
  parts <- unlist(strsplit(x, "/"))
  as.numeric(parts[1]) / as.numeric(parts[2])
})
ego_DKOvsWT_df$BgRatioNum <- sapply(ego_DKOvsWT_df$BgRatio, function(x) {
  parts <- unlist(strsplit(x, "/"))
  as.numeric(parts[1]) / as.numeric(parts[2])
})

#------------------------------------
# 4. Sort and select top 15 GO terms
#------------------------------------
ego_top_DKOvsWT <- ego_DKOvsWT@result %>%
  arrange(desc(Count)) %>%    # largest gene counts first
  slice(1:15)

# Wrap long GO term descriptions
ego_top_DKOvsWT$Description_wrapped <- str_wrap(ego_top_DKOvsWT$Description, width = 40)

# Set factor levels so largest dots appear at top
ego_top_DKOvsWT$Description_wrapped <- factor(
  ego_top_DKOvsWT$Description_wrapped,
  levels = rev(ego_top_DKOvsWT$Description_wrapped)
)

#------------------------------------
# 5. Visualize top 15 GO terms
#------------------------------------
ggplot(ego_top_DKOvsWT, aes(x = GeneRatio, y = Description_wrapped)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_continuous(low = "red", high = "blue", trans = "reverse") +
  scale_size(range = c(3, 8)) +
  theme_bw() +
  labs(
    title = "Top 15 enriched GO terms (DKO case vs WT case)",
    x = "Gene Ratio",
    y = NULL,
    color = "Adjusted p-value",
    size = "Gene count"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10, hjust = 1),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # vertical x-axis labels
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )

