#!/bin/bash
#SBATCH --job-name=CHIPseq
#SBATCH --nodes=2
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=6
#SBATCH --partition=2004
#SBATCH --mem=50gb
#SBATCH --output=pipeline_%j.log # Standard output and error log

##########################################################################################################################################################################
############################## peak calling sbatch script froma sample sheet #############################################################################
# This code will gather all the bam file names listed in the bam sample sheet from the input folder into an array to initiate pipeline runs for each one in the HPC
# The control bam is entered at the top of the sample sheet
# *****This script ideally should be initiated from the main pipeline but this requires further work**********
# The script creates a parent directory with the run ID, this is a date/time stamp unless specified as an option
# The script will create a folder for each bam file with the bam id as it's name within the parent folder
# Remember to change the SBATCH options above to configure for your run
# This script should only be run on the HPC using sbatch
# Options include:
#      threads = Will multi-thread any process to this number - Not currently used
#      input = Change the path of the input fastq files, default is ~/data
#      id = Change the name of the parent pipeline output folder, the default is a datestamp
# This only uses the default macs2 settings, options can be added as required
# If this is intitated from the main script, you would not reauire the while loop, maybe an if statement for how it is being used
# Author Steve Walsh May 2024
###########################################################################################################################################################################

#Uncomment below for debugging
#set -x

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

analysis_out_dir=${outdir}/${RUNID}
mkdir ${analysis_out_dir}/
mkdir ${analysis_out_dir}/MACS2

cp sample_sheet_macs2.csv ${analysis_out_dir}/sample_sheet_macs2_${RUNID}.csv
INPUT=${analysis_out_dir}/sample_sheet_macs2_${RUNID}.csv

#Get the control bam file from the first line of bam sample sheet
CONTROL_FILE=$(head -n 1 $INPUT)
sed -i '1d' $INPUT

while IFS= read -r LINE 
do

    # split line into array using , delimitator
    #echo ${var1}
    ARRAYLINE=(${LINE//,/ })
    BAM_FILE=${ARRAYLINE[0]}

    #Make directories for the peak call fileds as we go
    mkdir ${analysis_out_dir}/MACS2/${BAM_FILE}
    cd ${analysis_out_dir}/MACS2/${BAM_FILE}

    if [[ ${CONTROL_FILE} == "none" ]]; then
        srun --mem=10000MB --cpus-per-task=6 --ntasks=1 macs2 callpeak  -t ${bam_dir}/${BAM_FILE}.sorted.bam -f BAM -n ${BAM_FILE} &
    else
        srun --mem=10000MB --cpus-per-task=6 --ntasks=1 macs2 callpeak -c ${bam_dir}/${CONTROL_FILE}.sorted.bam -t ${bam_dir}/${BAM_FILE}.sorted.bam -f BAM -n ${BAM_FILE} &
    fi

done < ${INPUT}

#Wait for all pipelines to finish
wait
