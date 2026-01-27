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
sampleDists <- dist(t(assay(vsd)))
pheatmap(as.matrix(sampleDists))

#-----------------------------
# Checking groups and contrast
#-----------------------------

# check groups
design = ~ group
levels(dds$group)

# what pairwise comparisons are avaiable
resultsNames(dds)

# "Intercept"                                 
# "group_Blood_DKO_Control_vs_Blood_DKO_Case"  ok
# "group_Blood_WT_Case_vs_Blood_DKO_Case"      ok
# "group_Blood_WT_Control_vs_Blood_DKO_Case"  not ok 
#  I want:
# "group_Blood_WT_Control_vs_Blood_WT_Case"  

#----------
# Contrasts
#----------

# Blood_WT_Case vs Blood_DKO_Case
res_WT_case_vs_DKO_case <- results(
  dds,
  contrast = c("group", "Blood_WT_Case", "Blood_DKO_Case")
)

# Blood_DKO_Control vs Blood_DKO_Case
res_DKO_control_vs_DKO_case <- results(
  dds,
  contrast = c("group", "Blood_DKO_Control", "Blood_DKO_Case")
)

# Blood_WT_Control vs Blood_WT_Case
res_WT_control_vs_WT_case <- results(
  dds,
  contrast = c("group", "Blood_WT_Control", "Blood_WT_Case")
)

#--------------------
# Look at the results
#--------------------
#---------------------------------
# Blood_WT_Case vs Blood_DKO_Case
#---------------------------------

# View first few genes
head(res_WT_case_vs_DKO_case)

# Open in RStudio table viewer
View(res_WT_case_vs_DKO_case)

# Column names
colnames(res_WT_case_vs_DKO_case)

#-------------------------------
# Filter significant genes
# padj < 0.05 and |log2FC| > 1
#-------------------------------
res_sig <- res_WT_case_vs_DKO_case[
  !is.na(res_WT_case_vs_DKO_case$padj) &
    res_WT_case_vs_DKO_case$padj < 0.05 &
    abs(res_WT_case_vs_DKO_case$log2FoldChange) > 1,
]

# Number of significant genes
nrow(res_sig)

# Up- vs down-regulated genes
up_genes   <- sum(res_sig$log2FoldChange > 0)  # higher in WT_case
down_genes <- sum(res_sig$log2FoldChange < 0)  # higher in DKO_case

up_genes
down_genes

#-------------------------------
# Summary (all genes, for context)
#-------------------------------
summary(res_WT_case_vs_DKO_case, alpha = 0.05)

#-------------------------------
# MA plot
#-------------------------------
plotMA(
  res_WT_case_vs_DKO_case,
  xlab = "Mean normalized expression",
  ylab = "Log2 fold change",
  main = "WT case vs DKO case"
)

#-------------------------------
# Save significant genes to CSV
#-------------------------------
write.csv(
  as.data.frame(res_sig),
  file = "WT_case_vs_DKO_case_DESeq2_sig_genes.csv"
)

#-------------------------------------
# Blood_DKO_Control vs Blood_DKO_Case
#-------------------------------------

# View first few genes
head(res_DKO_control_vs_DKO_case)

# Open in RStudio table viewer
View(res_DKO_control_vs_DKO_case)

# Column names
colnames(res_DKO_control_vs_DKO_case)

#-------------------------------
# Filter significant genes
# padj < 0.05 and |log2FC| > 1
#-------------------------------
res_sig_DKO_control_vs_DKO_case <- res_DKO_control_vs_DKO_case[
  !is.na(res_DKO_control_vs_DKO_case$padj) &
    res_DKO_control_vs_DKO_case$padj < 0.05 &
    abs(res_DKO_control_vs_DKO_case$log2FoldChange) > 1,
]

# Number of significant genes
nrow(res_sig_DKO_control_vs_DKO_case)

# Up- vs down-regulated genes
up_genes_DKO   <- sum(res_sig_DKO_control_vs_DKO_case$log2FoldChange > 0)  # higher in DKO_Control
down_genes_DKO <- sum(res_sig_DKO_control_vs_DKO_case$log2FoldChange < 0)  # higher in DKO_Case

up_genes_DKO
down_genes_DKO

#-------------------------------
# Summary (all genes, for context)
#-------------------------------
summary(res_DKO_control_vs_DKO_case, alpha = 0.05)

#-------------------------------
# MA plot
#-------------------------------
plotMA(
  res_DKO_control_vs_DKO_case,
  xlab = "Mean normalized expression",
  ylab = "Log2 fold change",
  main = "DKO control vs DKO case"
)

#-------------------------------
# Save significant genes to CSV
#-------------------------------
write.csv(
  as.data.frame(res_sig_DKO_control_vs_DKO_case),
  file = "DKO_control_vs_DKO_case_DESeq2_sig_genes.csv"
)

#-------------------------------------
# Blood_WT_Control vs Blood_WT_Case
#-------------------------------------

# View first few genes
head(res_WT_control_vs_WT_case)

# Open in RStudio table viewer
View(res_WT_control_vs_WT_case)

# Column names
colnames(res_WT_control_vs_WT_case)

#-------------------------------
# Filter significant genes
# padj < 0.05 and |log2FC| > 1
#-------------------------------
res_sig_WT_control_vs_WT_case <- res_WT_control_vs_WT_case[
  !is.na(res_WT_control_vs_WT_case$padj) &
    res_WT_control_vs_WT_case$padj < 0.05 &
    abs(res_WT_control_vs_WT_case$log2FoldChange) > 1,
]

# Number of significant genes
nrow(res_sig_WT_control_vs_WT_case)

# Up- vs down-regulated genes
up_genes_WT   <- sum(res_sig_WT_control_vs_WT_case$log2FoldChange > 0)  # higher in WT_Control
down_genes_WT <- sum(res_sig_WT_control_vs_WT_case$log2FoldChange < 0)  # higher in WT_Case

up_genes_WT
down_genes_WT

#-------------------------------
# Summary (all genes, for context)
#-------------------------------
summary(res_WT_control_vs_WT_case, alpha = 0.05)

#-------------------------------
# MA plot
#-------------------------------
plotMA(
  res_WT_control_vs_WT_case,
  xlab = "Mean normalized expression",
  ylab = "Log2 fold change",
  main = "WT control vs WT case"
)

#-------------------------------
# Save significant genes to CSV
#-------------------------------
write.csv(
  as.data.frame(res_sig_WT_control_vs_WT_case),
  file = "WT_control_vs_WT_case_DESeq2_sig_genes.csv"
)

#----------------------------
# Functoin for Venn diagrams
#----------------------------
install.packages("VennDiagram")
install.packages("gridExtra")  # usually already installed, but safe


# Significant DE genes
genes_WTcase_vs_DKOcase <- rownames(res_WT_case_vs_DKO_case[
  !is.na(res_WT_case_vs_DKO_case$padj) &
    res_WT_case_vs_DKO_case$padj < 0.05 &
    abs(res_WT_case_vs_DKO_case$log2FoldChange) > 1,
])

genes_WTctrl_vs_WTcase <- rownames(res_WT_control_vs_WT_case[
  !is.na(res_WT_control_vs_WT_case$padj) &
    res_WT_control_vs_WT_case$padj < 0.05 &
    abs(res_WT_control_vs_WT_case$log2FoldChange) > 1,
])

genes_DKOctrl_vs_DKOcase <- rownames(res_DKO_control_vs_DKO_case[
  !is.na(res_DKO_control_vs_DKO_case$padj) &
    res_DKO_control_vs_DKO_case$padj < 0.05 &
    abs(res_DKO_control_vs_DKO_case$log2FoldChange) > 1,
])

library(VennDiagram)

venn.plot <- venn.diagram(
  x = list(
    WT_ctrl_vs_WT_case = genes_WTctrl_vs_WTcase,
    DKO_ctrl_vs_DKO_case = genes_DKOctrl_vs_DKOcase,
    WT_case_vs_DKO_case = genes_WTcase_vs_DKOcase
  ),
  filename = NULL,
  fill = c("skyblue", "pink", "lightgreen"),
  alpha = 0.6,
  cex = 1.2,
  cat.cex = 1.2
)

grid::grid.draw(venn.plot)

#-------------------------------
# UP-regulated genes
# padj < 0.05 & log2FC > 1
#-------------------------------
venn_up <- venn.diagram(
  x = list(
    "WT control vs. WT case"   = up_WTctrl_vs_WTcase,
    "DKO control vs. DKO case" = up_DKOctrl_vs_DKOcase,
    "WT case vs DKO case"  = up_WTcase_vs_DKOcase
  ),
  filename = NULL,
  fill = c("skyblue", "pink", "lightgreen"),
  alpha = 0.6,
  cex = 3,
  cat.cex = 1.8,
  main = "Up-regulated genes",
  main.cex = 2,
  cat.pos = c(-20, 20, 180),  # angles of labels for each set
  cat.dist = c(0.05, 0.05, 0.05) # distance from circle
)

grid::grid.newpage()
grid::grid.draw(venn_up)

#-------------------------------
# DOWN-regulated genes
# padj < 0.05 & log2FC > 1
#-------------------------------

down_WTctrl_vs_WTcase <- rownames(res_WT_control_vs_WT_case[
  !is.na(res_WT_control_vs_WT_case$padj) &
    res_WT_control_vs_WT_case$padj < 0.05 &
    res_WT_control_vs_WT_case$log2FoldChange < -1,
])

down_DKOctrl_vs_DKOcase <- rownames(res_DKO_control_vs_DKO_case[
  !is.na(res_DKO_control_vs_DKO_case$padj) &
    res_DKO_control_vs_DKO_case$padj < 0.05 &
    res_DKO_control_vs_DKO_case$log2FoldChange < -1,
])

down_WTcase_vs_DKOcase <- rownames(res_WT_case_vs_DKO_case[
  !is.na(res_WT_case_vs_DKO_case$padj) &
    res_WT_case_vs_DKO_case$padj < 0.05 &
    res_WT_case_vs_DKO_case$log2FoldChange < -1,
])


venn_down <- venn.diagram(
  x = list(
    "WT control vs. WT case"   = down_WTctrl_vs_WTcase,
    "DKO control vs. DKO case" = down_DKOctrl_vs_DKOcase,
    "WT case vs. DKO case"  = down_WTcase_vs_DKOcase
  ),
  filename = NULL,
  fill = c("skyblue", "pink", "lightgreen"),
  alpha = 0.6,
  cex = 3,
  cat.cex = 1.8,
  main = "Down-regulated genes",
  main.cex = 2,
  cat.pos = c(-20, 20, 180),  # angles of labels for each set
  cat.dist = c(0.05, 0.05, 0.05) # distance from circle
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

# sanity check
genes_of_interest_ensembl %in% rownames(res_WT_control_vs_WT_case)
# True gene will be labelled
# False geneis not present 

# check if DESeq row names are Ensembl IDs
head(rownames(res_WT_control_vs_WT_case))


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

# 1. WT Control vs WT Case
volcano_WT <- plot_volcano_ensembl(
  res = res_WT_control_vs_WT_case,
  title = "WT control vs WT case",
  genes_of_interest = genes_of_interest_ensembl,
  label_df = label_map
)

volcano_WT

# 2. DKO control vs DKO case

volcano_DKO <- plot_volcano_ensembl(
  res_DKO_control_vs_DKO_case,
  title = "DKO control vs DKO case",
  genes_of_interest = genes_of_interest_ensembl,
  label_df = label_map
)
volcano_DKO

# 3. WT case vs DKO case

volcano_WTvsDKO <- plot_volcano_ensembl(
  res_WT_case_vs_DKO_case,
  title = "WT case vs DKO case",
  genes_of_interest = genes_of_interest_ensembl,
  label_df = label_map
)
volcano_WTvsDKO

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

#====================================================
# Overrepresentation Analysis: WT case vs WT control
#====================================================

#----------------------------
# 2. Prepare gene lists
#----------------------------

# DE genes from previous step (filtered for padj < 0.05 & |log2FC| > 1)
DE_genes <- rownames(res_sig_WT_control_vs_WT_case)
all_genes <- rownames(res_WT_control_vs_WT_case)

# Remove any NA IDs
DE_genes <- DE_genes[!is.na(DE_genes)]
all_genes <- all_genes[!is.na(all_genes)]

# Optional checks
cat("Number of DE genes:", length(DE_genes), "\n")
cat("Number of all genes:", length(all_genes), "\n")
head(DE_genes)
head(all_genes)

#------------------------------------
# 3. Run overrepresentation analysis
#------------------------------------
ego <- enrichGO(
  gene          = DE_genes,         # DE genes
  universe      = all_genes,        # background genes
  OrgDb         = org.Mm.eg.db,     # mouse annotation database
  ont           = "ALL",             # GO category: BP, MF, CC, or ALL
  keyType       = "ENSEMBL",        # our IDs are Ensembl
  pAdjustMethod = "BH",             # multiple testing correction
  pvalueCutoff  = 0.05,             # significance threshold
  qvalueCutoff  = 0.05,             # FDR threshold
  readable      = TRUE              # converts Ensembl IDs to gene symbols
)

#----------------------------------
# 4. Convert results to data frame
#----------------------------------
ego_df <- as.data.frame(ego)

#------------------------------------------------
# 5. Calculate GeneRatio and BgRatio numerically
#------------------------------------------------
# GeneRatio = DE genes in term / total DE genes
ego_df$GeneRatioNum <- sapply(ego_df$GeneRatio, function(x) {
  parts <- unlist(strsplit(x, "/"))
  as.numeric(parts[1]) / as.numeric(parts[2])
})

# BgRatio = all genes in term / total background genes
ego_df$BgRatioNum <- sapply(ego_df$BgRatio, function(x) {
  parts <- unlist(strsplit(x, "/"))
  as.numeric(parts[1]) / as.numeric(parts[2])
})

#----------------------------------------------
# 6. Sort by adjusted p-value & extract top 10
#----------------------------------------------
top10_GO <- ego_df[order(ego_df$p.adjust), ][1:10,
                                             c("ID", "Description", "GeneRatio", "GeneRatioNum",
                                               "BgRatio", "BgRatioNum", "p.adjust")]
top10_GO

#----------------------------
# 7. Visualize results
#----------------------------
# Dotplot
dotplot(ego, showCategory=10, title="Top 10 enriched gene ontology terms for WT control vs WT case") +
  theme(axis.text.y=element_text(size=10))

#=======================================================
# Overrepresentation Analysis: DKO control vs DKO case
#=======================================================

#----------------------------
# 1. Prepare gene lists
#----------------------------
# DE genes for DKO control vs DKO case comparison
DE_genes_DKO <- rownames(res_sig_DKO_control_vs_DKO_case)
all_genes_DKO <- rownames(res_DKO_control_vs_DKO_case)

# Remove any NA IDs
DE_genes_DKO <- DE_genes_DKO[!is.na(DE_genes_DKO)]
all_genes_DKO <- all_genes_DKO[!is.na(all_genes_DKO)]

# Optional checks
cat("Number of DE genes (DKO):", length(DE_genes_DKO), "\n")
cat("Number of all genes (DKO):", length(all_genes_DKO), "\n")
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
# 4. Sort and report top 10 GO terms
#------------------------------------
top10_GO_DKO <- ego_DKO_df[order(ego_DKO_df$p.adjust), ][1:10,
                                                         c("ID", "Description", "GeneRatio", "GeneRatioNum",
                                                           "BgRatio", "BgRatioNum", "p.adjust")]
top10_GO_DKO

#----------------------------
# 5. Visualize results
#----------------------------
dotplot(ego_DKO, showCategory=10, title="Top 10 enriched gene ontology terms for DKO control vs DKO case") +
  theme(axis.text.y=element_text(size=10))

#==================================================
# Overrepresentation Analysis: WT case vs DKO case
#==================================================

#-----------------------------------------------
# 1. Prepare gene lists for WT case vs DKO case
#-----------------------------------------------
# DE genes for WT case vs DKO case comparison
DE_genes_case <- rownames(res_sig_WT_case_vs_DKO_case)
all_genes_case <- rownames(res_WT_case_vs_DKO_case)

# Remove any NA IDs
DE_genes_case <- DE_genes_case[!is.na(DE_genes_case)]
all_genes_case <- all_genes_case[!is.na(all_genes_case)]

# Optional checks
cat("Number of DE genes (WT vs DKO case):", length(DE_genes_case), "\n")
cat("Number of all genes (WT vs DKO case):", length(all_genes_case), "\n")
head(DE_genes_case)
head(all_genes_case)

#------------------------------------
# 2. Run overrepresentation analysis
#------------------------------------
ego_case <- enrichGO(
  gene          = DE_genes_case,      # DE genes
  universe      = all_genes_case,     # background genes
  OrgDb         = org.Mm.eg.db,       # mouse annotation database
  ont           = "ALL",               # GO category
  keyType       = "ENSEMBL",          # our IDs are Ensembl
  pAdjustMethod = "BH",               # multiple testing correction
  pvalueCutoff  = 0.05,               # significance threshold
  qvalueCutoff  = 0.05,               # FDR threshold
  readable      = TRUE                # converts Ensembl IDs to gene symbols
)

#----------------------------------
# 3. Convert results to data frame
#----------------------------------
ego_case_df <- as.data.frame(ego_case)

# Calculate numeric GeneRatio and BgRatio
ego_case_df$GeneRatioNum <- sapply(ego_case_df$GeneRatio, function(x) {
  parts <- unlist(strsplit(x, "/"))
  as.numeric(parts[1]) / as.numeric(parts[2])
})
ego_case_df$BgRatioNum <- sapply(ego_case_df$BgRatio, function(x) {
  parts <- unlist(strsplit(x, "/"))
  as.numeric(parts[1]) / as.numeric(parts[2])
})

#------------------------------------
# 4. Sort and report top 10 GO terms
#------------------------------------
top10_GO_case <- ego_case_df[order(ego_case_df$p.adjust), ][1:10,
                                                            c("ID", "Description", "GeneRatio", "GeneRatioNum",
                                                              "BgRatio", "BgRatioNum", "p.adjust")]
top10_GO_case

#----------------------------
# 5. Visualize results
#----------------------------
dotplot(ego_case, showCategory=10, title="Top 10 enriched gene ontology terms for WT case vs DKO case") +
  theme(axis.text.y=element_text(size=10))




