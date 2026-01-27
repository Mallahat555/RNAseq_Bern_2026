This is step by step guidance on how to analyse tables obtained from featureCounts in step 4.
You don't need to work on the cluster for this part. 

1. Install software on you local machine:

    1.1 R version 4.5.1 (2025-06-13 ucrt) -- "Great Square Root"
    1.2 Bioconductor version: release 3.22 (BiocManager 1.30.27) 
        DESeq2 version: 1.50.2 with packages:
            Packages and their versions:
                pheatmap_1.0.13
                ggplot2_4.0.1
                DESeq2_1.50.2
               

2. Make DESeq2-readable table:

    For comfort I made DESeq2 tables while still on cluster and downloaded the output file. But you can download the 
    featureCounts reads_per_gene.txt files, which are in /data/users/apiatkowska/RNAseq/Step_4_count_reads/reads_per_gene_per_sample
    directory and process in on you local machine.

    DESeq2-readable table contains Gene ID and read counts per sample for all samples.
    
    2.1 Use this script to generate the tables per sample and a combined table for all samples:
        make_deseq2_counts.sh
    2.2 change the file permissions with "chmod": make the file executable "+x":
        in terminal put: chmod +x make_deseq2_counts.sh 
    2.3 run script in the terminal: ./make_deseq2_counts.sh
        make sure you are in the right directory to execute the script
    2.4 inspect resulting counts_all_samples.txt and *_counts.txt files for errors
    2.5 you will work with counts_all_samples.txt for the next step:
        /data/users/apiatkowska/RNAseq/Step_5_data_analysis/tables_of_counts_for_DESeq2/counts_all_samples.txt

3.  Work in R with DESeq2 installed on your local machine.
    3.1 Download the DESeq2-readable tables to your local machine. 
    3.2 See R script for details on how the analysis was performed:
        DESeq2_analysis.R
    
    NOTE! I defined a path with counts_all_samples.txt on my local machine and set it up to a variable counts_file
    The entire R code relies on that variable in further analysis.
    To repeat my analysis make sure you defined path specific to your machine and set up the same variable name. 
    If variable names differs, you must chabge the name in the code which follows. It best to keep my variable name.
    counts_file <- "C:/Users/agnie/Desktop/MSc Bioinformatics/MSc Bioinformatics Courses/RNA-seq/Project_RNA_seq/Step_5_data_analysis/tables_of_counts_for_DESeq2/counts_all_samples.txt"
    
    

