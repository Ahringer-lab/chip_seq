#!/bin/bash
#SBATCH --job-name=RNASeq  
#SBATCH --nodes=6
#SBATCH --ntasks=30
#SBATCH --cpus-per-task=20
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
#      jobs = Changes the maximum number of jobs that are submitted to the cluster at once, this option should be used in conjunction with the slurm options above
#      star_index = The location of the STAR index
#      kallisto_index = The location of the Kallisto index
# Author Steve Walsh May 2024
###########################################################################################################################################################################

set -x

#Set the defaults
outdir=~/out
fastq_dir=~/data/
BOWTIE_INDEX=~/references/built_genomes/star/c.elegans.latest
CHROM_SIZES=/mnt/home3/ahringer/index_files/genomes/c_elegans.PRJNA13758.WS285.genomic.chrom.sizes
THREADS=1
RUNID="PipelineRun-$(date '+%Y-%m-%d-%R')"
MERGEID=merged
WD="$(pwd)"
JOBS=1

# Function to handle incorrect arguments
function exit_with_bad_args {
    echo "Usage: bash lane_merger.bash optional args: --threads <number of threads> --input <input path> --id <Run ID>  --mergeID <merge ID> --jobs <Number of pairs of fastqs to run> --bowtie_index"
    echo "Invalid arguments provided" >&2
    exit # this stops the terminal closing when run as source
}

#Set the possible input options
options=$(getopt -o '' -l threads: -l input: -l id: -l mergeID -l jobs: -l bowtie_index -- "$@") || exit_with_bad_args

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
        --jobs)
            shift
            JOBS="$1"
            ;;
        --bowtie_index)
            shift
            BOWTIE_INDEX="$1"
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
echo "$analysis_out_dir"

#Set up stats folder
mkdir ${analysis_out_dir}/stats

MERGEID=_merged

cd ${WD}

#Set sample sheet name and a counter for number of jobs to be sent to the hpc at once
INPUT="sample_sheet.csv"
sed -i '/^[[:space:]]*$/d' $INPUT
sed -i 's/\r$//' $INPUT

COUNTER=0

while IFS= read -r LINE 
do

    # split line into array using tab delimitator - 0: FASTQ: FASTQ FILE SAMPLE_NAME: Name of sample (If changing from fastq file name)
    #echo ${var1}
    ARRAYLINE=(${LINE//,/ })
    FASTQ=${ARRAYLINE[0]}
    SAMPLE_NAME=${ARRAYLINE[1]}

#Loops through the fastq names, make directories for each output, ${base} holds the sample id (TODO Chane $base to something else)
    echo "Fastq file being analysed"
    echo "${FASTQ}${MERGEID}_R1_001.fastq.gz"
    echo "${FASTQ}${MERGEID}_R2_001.fastq.gz"
    echo "Sample ID being used"
    echo ${SAMPLE_NAME}

    srun --mem=10000MB --cpus-per-task=6 --ntasks=1 ./chip_seq_pipeline.bash --fastqid ${FASTQ} --sample_id ${SAMPLE_NAME} --threads ${THREADS} --input ${fastq_dir} --id ${RUNID} --mergeID ${MERGEID} --bowtie_index ${BOWTIE_INDEX} &

    #Jobs are submitted to the hpc up to the upper job limit.
    COUNTER=$(( COUNTER + 1 ))
    echo $COUNTER
    if [ "$COUNTER" -ge "$JOBS" ]; then
        echo "Maximum number of pipelines are running ($JOBS), waiting for them to finish"
        wait
        unset COUNTER
        echo "Running the next group of pipelines now"
        COUNTER=0
    fi

done < ${INPUT}

#Wait for all pipelines to finish
wait

#Carry out multiqc across all samples
cd ${analysis_out_dir}
multiqc .

Make the summary stats file
cd ${analysis_out_dir}/stats
echo \#Run ID\: ${RUNID} > summary_stats.txt
echo -n Sample_name, >> summary_stats.txt
awk 'FNR==3{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Forward_fastq, >> summary_stats.txt
awk 'FNR==4{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Reverse_fastq, >> summary_stats.txt
awk 'FNR==5{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Aligned_reads, >> summary_stats.txt
awk 'FNR==6{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Q30_Aligned_reads, >> summary_stats.txt
awk 'FNR==7{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Q10_Aligned_reads, >> summary_stats.txt
awk 'FNR==8{printf $0 >> "summary_stats.txt"}' *
sed 's/.$//' summary_stats.txt >> summary_stats.csv
rm summary_stats.txt
