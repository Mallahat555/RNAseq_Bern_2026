This is step by step description on how to map reads.

Species: Mus muluscus
Genome version: GRCm39
Gene annotation release: 115
Ensembly file: Mus_musculus.GRCm39.115.gtf.gz
HISAT2 version 2.2.1 
Samtools version 1.20

1. Preparation
    1.1 Prepare appropriate project tree:
        My tree is in /data/users/apiatkowska/RNAseq/   directory:

        Step_3_map_reads/
            reference_genome_and_annotations/
            hisat2_indexing/
            mapped_reads/
                --> this part of the tree is created by run_map_reads_for...slurm code
                sam/
                bam/
                bam_sorted/
                logs/ 
            mapped_reads_small_regions_IGV/
                bam_sorted_small_reads/
                logs/

        useful command lines:
        mkdir -- make directory
        ls -- to check what is in the folder
        cd ./ -- go to the appropriate directory:
        cd ./data/users/apiatkowska/RNAseq/Step_3_map_reads/reference_genome_and_annotations

    1.2 Containers directories:
        hisat2: /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif

2. Download reference genome and annotations from Ensembly:
    - to see the directory content of Ensembly:
    wget -qO- ftp://ftp.ensembl.org/pub/current_gtf/mus_musculus/ | less

    2.1 Download reference genome FASTA:
        wget ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

        unzip fasta file:
        gunzip -k /data/users/apiatkowska/RNAseq/Step_3_map_reads/reference_genome_and_annotations/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
        -k keeps the original zipped file instead of deleting it

    2.2 Download annotation GTF:
        wget ftp://ftp.ensembl.org/pub/current_gtf/mus_musculus/Mus_musculus.GRCm39.115.gtf.gz

        unzip GTF file:
        gunzip -k /data/users/apiatkowska/RNAseq/Step_3_map_reads/reference_genome_and_annotations/Mus_musculus.GRCm39.115.gtf.gz

        ✔ Full Ensembl gene annotation
        ✔ Includes all chromosomes + scaffolds
        ✔ Standard choice for RNA-seq gene counting
        ✔ Matches genome assembly (GRCm39)

    2.3 Download checksums:
        wget ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/CHECKSUMS
        wget ftp://ftp.ensembl.org/pub/current_gtf/mus_musculus/CHECKSUMS

    2.4 Verify if the files are in your directory:
        ls 

    2.5 Verify checksums FASTA and GTF integrity:
        - Required step to know if both files GFT and FASTA are not corrupted.

        I verified GTF and obtained this results:
        sum Mus_musculus.GRCm39.115.gtf.gz
        Output: 
        48874 40705
        Compare to GTF file: 
        48874 40705 /hps/nobackup/flicek/ensembl/production/release_dumps/release-115/ftp_dumps/vertebrates/gtf/mus_musculus/Mus_musculus.GRCm39.115.gtf.gz Mus_musculus.GRCm39.115.gtf.gz

        I verified FASTA and obtained this results:
        sum Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
        Output: 16996 787519
        Compare to FASTA file: 
        16996 787519 /hps/nobackup/flicek/ensembl/production/release_dumps/release-115/ftp_dumps/vertebrates/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

        For both FASTA and GTF:
        ✔ Checksum value matches
        ✔ File size matches
        ✔ Filename matches
        ✔ Assembly matches (GRCm39)
        ✔ Ensembl release matches (115) 
        ✔ No corruption

3. Build HISAT2 index (once).
    - uses only FASTA (not the GTF), remember to unzip fasta file:
        gunzip -k /data/users/apiatkowska/RNAseq/Step_3_map_reads/reference_genome_and_annotations/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
    - change directory to cd ./data/users/apiatkowska/RNAseq/Step_3_map_reads/histat2_indexing

    3.1 Develop code in .slurm file:
        - see histat2_index.slurm (I developed my code there)
    3.2 Submitt the job:
        sbatch histat2_index.slurm
    3.3 Verify the index files
        ls /data/users/apiatkowska/RNAseq/Step_3_map_reads/hisat2_indexing/index/

4. Map reads to the reference genome.
    The code performs following tasks:

    Maps reads with HISTAT2 to the reference genome.
    The correct strandedness setting for this library prep protocol is RF.
    SAM files are converted to BAM files. 
    BAM files are sorted by genomic coordinates. 
    BAM files are indexed.

    4.1 The code to run just one sample and debug on demand.
        sbatch run_map_reads_for_one_file.slurm

    4.2 Modify the code to run on all samples.
        Here I used multiple arrays 0-14 instead of forward loop. It is faster, about 30 min run time.
        sbatch run_map_reads_for_all_samples.slurm

        useful command to check how samples are assigned per array:
        grep "Processing sample" /data/users/apiatkowska/RNAseq/Step_3_map_reads/mapped_reads/logs/*.out
    
5. Inspect BAM with IGV (optional)

    5.1 Extract a small genomic region from bam_sorted files for each biological sample, to reduce their size:
        selected region REGION="chr1:150000000-150100000"
        sbatch run_extract_small_regions_slurm
        and download them to your local machine

    5.2 Install IGV from https://igv.org/doc/desktop/#DownloadPage/

    5.3 In IGV:
        5.3.1 Load the correct reference genome:
        Download the reference genome to your machine, needs to be accessible locally, and unzipped. 
        In IGV top menue → Genomes → Load Genome from file → select where you stored the unzippes .fa file

        5.3.2 Load selected BAM file into IGV:
        File → Load from File
        Select:
        .bam IGV automatically looks for .bai file 

        5.3.3 Load annotation GTF:
        File → Load from File
        Mus_musculus.GRCm39.115.gtf
        This helps verify:
        Reads align to known genes
        Strandness looks correct (important for RF libraries)

    5.4 Inspect alignments:

        ✔ Reads pile up on exons
        ✔ Spliced reads (gaps = introns)
        ✔ Paired-end reads facing inward
        ✔ Consistent coverage
        ✔ No massive mismapping

        Zoom in/out to inspect exon–intron structure.


