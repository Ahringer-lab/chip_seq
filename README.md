## About

This repository contains a selection of shell scripts to carry out the primary analysis of CHIP or ATAC-seq data.

The main pipeline sbatch_main_CHIPseq_pipeline.bash will read in fastq files, carry out qc on the files, align the files to the c.elegans genome and then convert the BAM files to bigwigs.

The peakcalling pipeline can be used to call peaks using MACS2 on BAM files

### Install

To use these pipelines first clone the repository into your home directory on the cluster as follows:

git clone https://github.com/Ahringer-lab/bulk_rna_seq.git

Then create two directories in your home directory:

* data

* out

Finally clone the conda environment as below:

conda env create -f environment.yml

### Running the pipeline main pipeline

First switch to the correct conda environment:

conda activate CHIP-Seq

Drop all files into the data folder, by default the pipeline will take files with the following format:

<my_fastq>_merged_R1/2_001.fastq.gz

This will usually be files that have been merged with the file merger script, if you want to use a differnt file name format the mergeID flag should be used (see below)

From within the directory containing the scripts first edit the sample_sheet.csv file, you should list the files you want to analyse in the first column without _merged_R1/2_001.fastq.gz in the file name and 
in the second column list any updates you want to make to the output file names. 

Once the sample_sheet.csv file is ready fun the pipeline with the following command:

sbatch sbatch_main_CHIPseq_pipeline.bash

*Note: You may have to change the Slurm cluster settings at the top of sbatch_pipeline.bash depending on how many files you are analysing and how busy the cluster is.

The following flags can be used to change the behaviour of the pipeline:

--threads : To run the individual programs within the pipeline by this many threads (The input file are run in tandem by default)

--input : Change the input folder from the default

--id : Change the ouput folder name from the default of a date/time stamp

--mergeID : Change the default input from <my_fastq>_merged_R1/2_001.fastq.gz. It will change the _merged_ part in the middle.

--bwa_index : Change the dedault location of the bwa index from /mnt/home3/ahringer/index_files/built_indexes/bwa/c_elegans.PRJNA13758.WS285/c_elegans.PRJNA13758.WS285.genomic.fa

N.B The defaul index is built from the c_elegans.PRJNA13758.WS285.genomic.fa file from Wormbase

### Output

All files will be ouput to the out directory in your home directory, there will be a folder for each sample that contains the following folders:

* fastq  
* fastq_screen  
* bwa  
* trim_galore
* bw

N.B default settings are used for all programs and bigwigs are created with a bin size of 10 and normalised by counts per million (CPM), the reads are not shifted or extended.

The fastq folder contains a copy of the fastq files used (as a point of reference to check back) and the other folders contain the output of the respective programs.

### Running the peakcalling pipeline

The pipeline is run in the same way as the main pipeline except you put the bam files in the input folder rather than the fastq files and a different sample sheet is used. For this pipeline the sample sheet called sample_sheet_macs2.csv is used instead, it is set up in exactly the same way as the sample sheet for the main pipleine, but in the first row put the name of the control bam file if using one, if not put 'none' instead.

The flags for this pipeline include:

--threads : To run the individual programs within the pipeline by this many threads (The input file are run in tandem by default)

--input : Change the input folder from the default

--id : Change the ouput folder name from the default of a date/time stamp

### Additional information

There are several additional 'mini' pipelines in the repository including:

* bw.bash : This converts bam files to bw files, it works exactly like the main pipeline but with bams as the input and used the sample_sheet_bw.csv as its sample sheet.
* control.bash : This file aligns fastq files using bwa, this script had been intended to be used in the main pipeline if peak calling was integrated into it
* stats.bash : This is used to summarise all the individual stats file generated for each sample

