#!/bin/bash
#SBATCH --job-name=CHIPseq
#SBATCH --nodes=2
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=6
#SBATCH --partition=2004
#SBATCH --mem=50gb
#SBATCH --output=pipeline_%j.log # Standard output and error log

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
mkdir ${analysis_out_dir}/BW

while IFS= read -r LINE 
do

    # split line into array using , delimitator
    #echo ${var1}
    ARRAYLINE=(${LINE//,/ })
    BAM_FILE=${ARRAYLINE[0]}

    #Make directories for the peak call fileds as we go
    mkdir ${analysis_out_dir}/BW/${BAM_FILE}

    bamCoverage -p ${THREADS} --normalizeUsing CPM -bs 50 -b ${bam_dir}/${BAM_FILE}.sorted.bam -o ${analysis_out_dir}/bw/${BAM_FILE}.bw

done < ${INPUT}