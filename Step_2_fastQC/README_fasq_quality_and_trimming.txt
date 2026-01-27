This is step by step description on how to check the quality of raw fastq reads, trimm them and check the quality again.

fastQC version: 0.12.1
MultiQC version: 1.19
fastP version: 0.24.1

1. Preparation:

   1.1 Prepare appropriate project tree.
   My tree is in /data/users/apiatkowska/RNAseq/   directory:

    Step_2_fastQC/
        fastQC_results/
            error/
            output/
        multiQC_results/
    Step_2.2_trimmed_fastQ/
        individual_fastP_quality_control_reports/
            error/
            output/
        multiQC_fastp_results/
            error/
            output/
        trimmed_fastQ/

    useful commands: 

    mkdir -- to make a folder
    ls -- to check what is in the folder
    cd ./ -- go to the appropriate directory:

    1.2 Useful Directories:

        For raw fastQ files: 
        /data/courses/rnaseq_course/toxoplasma_de/reads_Blood

        For containers:
        fastQC: /containers/apptainer/fastqc-0.12.1.sif
        MultiQC: /containers/apptainer/multiqc-1.19.sif
        FastP: /containers/apptainer/fastp_0.24.1.sif

2. Run quality control on raw fastQ files:

    I developed myh codes in .sh files, as job is small.
    
    For indivudual raw fastQ files:
    Go to the appropriate directory 
    /data/users/apiatkowska/RNAseq/Step_2_fastQC 
    and submit the job from the terminal
    sbatch run_raw_fastQ_quality_control.sh

    For multiple fastQ files:
    In the same directory, submitt job:
    sbatch run_raw_fastQ_MultiQC.sh

3. Trim fastQ files:

    Trimming NEEDS TO BE PAIR ENDED!!!

    Go to the appropriate directory 
    /data/users/apiatkowska/RNAseq/Step_2.2_trimmed_fastQ
    and submit the job from the terminal
    sbatch run_fastQ_trimming_pair_ended.sh

4. Quality contor on trimmed fastQ files:

    Go to the appropriate directory 
    /data/users/apiatkowska/RNAseq/Step_2.2_trimmed_fastQ
    and submit the job from the terminal

    For individual fastP files:
    sbatch run_fastP_quality_control.sh

    For multiple trimmed fastQ files:
    sbatch run_trimmed_fastQ_MultiQC_repport.sh

You are done with this step.



            
        
     