#!/bin/bash
#SBATCH --job-name=CHIPSeq  
#SBATCH --nodes=6
#SBATCH --ntasks=60
#SBATCH --cpus-per-task=12
#SBATCH --mem=100gb
#SBATCH --partition=2004
#SBATCH --output=pipeline_%j.log # Standard output and error log

##########################################################################################################################################################################
############################## bulk rna-seq sbatch job submission script with a sample sheet #############################################################################
# This code will gather all the fastq file names listed in the sample sheet from the input folder into an array to initiate pipeline runs for each one in the HPC
# Using the sample sheet the output file and folders can be named differently to the fastq file names, this is useful if the fastq files have been named non-intuitively
# The pipeline assumes fastq file have been merged using the bash script on github
# The fastq id is the first part of the standard file name afer merging e.g. for JAtab71-F-1_merged_R1_001.fastq.gz the id is JAtab71-F-1
# The sample ID can be set to whatever the user wants
# The script creates a parent directory with the run ID, this is a date/time stamp unless specified as an option
# The script will create a folder for each pair of fastq files with the fastq id as it's name within the parent folder
# Remember to change the SBATCH options above to configure for your run, ntasks should be the number of fastq pairs
# This script should only be run on the HPC using sbatch
# Options include:
#      threads = Will multi-thread any process to this number
#      input = Change the path of the input fastq files, default is ~/data
#      id = Change the name of the parent pipeline output folder, the default is a datestamp
#      mergeID = If the file names have been merged differently the input can be changed here 'fastqid_<Add the flag here>_R1/R2_001.fastq.gz'
# Currently this script will not proceed beyond the alignment step, I have added some groundwork for a full pipeline (Commented out)
# Peak calling is run in another script but could be initiated here at some point.
# The option to merge bam file could be added in and handled via the sample sheet
# Author Steve Walsh May 2024
###########################################################################################################################################################################

#Uncomment below for debugging
#set -x

#Set the defaults
outdir=~/out
fastq_dir=~/data/
BWA_INDEX=/mnt/home3/ahringer/index_files/built_indexes/bwa/c_elegans.PRJNA13758.WS285/c_elegans.PRJNA13758.WS285.genomic.fa
CHROM_SIZES=/mnt/home3/ahringer/index_files/genomes/c_elegans.PRJNA13758.WS285.genomic.chrom.sizes
THREADS=1
RUNID="PipelineRun-$(date '+%Y-%m-%d-%R')"
MERGEID=merged
WD="$(pwd)"
SAMPLE_NAME=null

# Function to handle incorrect arguments
function exit_with_bad_args {
    echo "Usage: bash chip_seq_alignment.bash optional args: --threads <number of threads> --input <input path> --id <Run ID> --mergeID <merge ID> --bwa_index"
    echo "Invalid arguments provided" >&2
    exit # this stops the terminal closing when run as source
}

#Set the possible input options
options=$(getopt -o '' -l threads: -l input: -l id: -l mergeID -l bwa_index: -- "$@") || exit_with_bad_args

#Get the inputs
eval set -- "$options"
while true; do
    case "$1" in
        --threads)
            shift
            THREADS="$1"
            ;;
        --input)
            shift
            fastq_dir="$1"
            ;;
        --id)
            shift
            RUNID="$1"
            ;;
        --mergeID)
            shift
            MERGEID="$1"
            ;;
        --bwa_index)
            shift
            BWA_INDEX="$1"
            ;;
         --)
            shift
            break
            ;;
    esac
    shift
done

#Set and create the ouput directory based on Run ID (Date/time if not set)
analysis_out_dir=${outdir}/${RUNID}
mkdir $analysis_out_dir

#Set sample sheet name and a counter for number of jobs to be sent to the hpc at once
INPUT="sample_sheet.csv"

#****The code below reads in an updated sample sheet with the control bam or fastq files identified in the first two lines and aligns the control bam files****
#This would be the first step in fully automating this process
#The pipeline should either take fastq file or a bam for the control

#Get the control file and determine if it is in bam or fastq format
#sed -i '/^[[:space:]]*$/d' $INPUT
#sed -i 's/\r$//' $INPUT
#CONTROL_TYPE=$(head -n 1 $INPUT)
#sed -i '1d' $INPUT
#CONTROL_FILE=$(head -n 1 $INPUT)
#sed -i '1d' $INPUT

#Check the sample sheet contains the control file
#if [ ${CONTROL_TYPE} != "bam" ] && [ ${CONTROL_TYPE} != "fastq" ]; then
#    echo "Sample sheet type not found or incorrect"
#    exit 1
#elif [ ${CONTROL_TYPE} == "null" ]; then
#    echo "Sample sheet type not found or incorrect"
##    exit 1
#elif [[ ${CONTROL_FILE} == 'null' ]]; then
#    echo "Sample sheet not found or incorrect"
#    exit 1
#fi

#Generate the control bam file for peak calling if fastq files are given as an input
#if [ ${CONTROL_TYPE} == "fastq" ]; then
#    srun --nodes=1 --mem=15000MB --cpus-per-task=6 --ntasks=1 ./control.bash --fastqid ${CONTROL_FILE} --threads ${THREADS} --input ${analysis_out_dir} --id ${RUNID}
#   CONTROL_FILE=${analysis_out_dir}/control/${CONTROL_FILE}.bam
#done

#From this point the script carries out quality control on the fastq files and aligns

#Set up stats folder
mkdir ${analysis_out_dir}/stats

MERGEID=_merged

cd ${WD}

cp ${INPUT} ${analysis_out_dir}/Sample_sheet_${RUNID}.log

while IFS= read -r LINE 
do

    # split line into array using tab delimitator - 0: FASTQ: FASTQ FILE SAMPLE_NAME: Name of sample (If changing from fastq file name)
    #echo ${var1}
    ARRAYLINE=(${LINE//,/ })
    FASTQ=${ARRAYLINE[0]}
    SAMPLE_NAME=${ARRAYLINE[1]}

    

#Loops through the fastq names, make directories for each output, ${base} holds the sample id (TODO Chane $base to something else)

    srun --mem=10000MB --cpus-per-task=6 --ntasks=1 ./chip_seq_alignment.bash --fastqid ${FASTQ} --sampleid ${SAMPLE_NAME} --threads ${THREADS} --input ${fastq_dir} --id ${RUNID} --mergeID ${MERGEID} --bwaindex ${BWA_INDEX} &

#Carry out peak calling
#Ideally the fastq file that have been aligned would be peak called here, this requires organising the control and treatment files correctly, see code above

# srun --nodes=1 --mem=1000MB --cpus-per-task=6 --ntasks=1 ./peakcalling.bash --input ${analysis_out_dir}/bam --id ${RUNID}

done < ${INPUT}

#Wait for all pipelines to finish
wait

#Carry out multiqc across all samples
cd ${analysis_out_dir}
multiqc .

#Make the summary stats file
cd ${WD}
srun --nodes=1 --mem=5000MB --cpus-per-task=6 --ntasks=1 ./stats.bash --input ${analysis_out_dir} --id ${RUNID}
