#!/bin/bash
#SBATCH --job-name=RNASeq  
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=12
#SBATCH --mem=50gb
#SBATCH --partition=2004
#SBATCH --output=pipeline_%j.log # Standard output and error log

##########################################################################################################################################################################
############################## bulk rna-seq sbatch job submission script with a sample sheet #############################################################################
# This code will gather all the bam file names listed in the sample sheet from the input folder into an array to initiate pipeline runs for each one in the HPC
# The control bam is entered at the top of the sample sheet
# *****This script ideally should be initiated from the main pipeline but this requires further work**********
# The script creates a parent directory with the run ID, this is a date/time stamp unless specified as an option
# The script will create a folder for each pair of fastq files with the fastq id as it's name within the parent folder
# Remember to change the SBATCH options above to configure for your run, ntasks should be the number of fastq pairs
# This script should only be run on the HPC using sbatch
# Options include:
#      threads = Will multi-thread any process to this number
#      input = Change the path of the input fastq files, default is ~/data
#      id = Change the name of the parent pipeline output folder, the default is a datestamp
#      mergeID = If the file names have been merged differently the input can be changed here 'fastqid_<Add the flag here>_R1/R2_001.fastq.gz'
# Author Steve Walsh May 2024
###########################################################################################################################################################################

set -x

#Set the defaults
outdir=~/out
bam_dir=~/data/
CHROM_SIZES=/mnt/home3/ahringer/index_files/genomes/c_elegans.PRJNA13758.WS285.genomic.chrom.sizes
THREADS=1
RUNID="PipelineRun-$(date '+%Y-%m-%d-%R')"
WD="$(pwd)"

# Function to handle incorrect arguments
function exit_with_bad_args {
    echo "Usage: bash lane_merger.bash optional args: --threads <number of threads> --input <input path> --id <Run ID>"
    echo "Invalid arguments provided" >&2
    exit # this stops the terminal closing when run as source
}

#Set the possible input options
options=$(getopt -o '' -l threads: -l input: -l id: -- "$@") || exit_with_bad_args

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
            bam_dir="$1"
            ;;
        --id)
            shift
            RUNID="$1"
            ;;
         --)
            shift
            break
            ;;
    esac
    shift
done

INPUT="sample_sheet_macs2.csv"

#Get the control bam file from the first line of bam sample sheet
CONTROL_FILE=$(head -n 1 $INPUT)
sed -i '1d' $INPUT

analysis_out_dir=${outdir}/${RUNID}

mkdir ${analysis_out_dir}/MACS2

while IFS= read -r LINE 
do

    # split line into array using , delimitator
    #echo ${var1}
    ARRAYLINE=(${LINE//,/ })
    BAM_FILE=${ARRAYLINE[0]}

    #Make directories for the peak call fileds as we go
    mkdir ${analysis_out_dir}/MACS2/${BAM_FILE}
    cd ${analysis_out_dir}/MACS2/${BAM_FILE}

    srun --mem=10000MB --cpus-per-task=6 --ntasks=1 macs2 callpeak -c ${bam_dir}/${CONTROL_FILE}.sorted.bam -t ${bam_dir}/${BAM_FILE}.sorted.bam -f BAM -n ${BAM_FILE} &

done < ${INPUT}

#Wait for all pipelines to finish
wait
