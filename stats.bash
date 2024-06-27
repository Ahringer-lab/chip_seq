#!/bin/bash

###############################################################################################################################################
############################## Stats summarizer ###############################################################################################
#This code will take the inputs of individual stats files from pipeline run and summarise them together
#The summary stats file will include the fastq read numbers and alignments stats for all individual stats files given as an input
#This script is run each time a batch pipeline run completes, but can be run any time with command line options
#Options include:
#   input = The directory containing the stats files to be summarized
#   id = The pipeline ID if summarising stats files from a previous pipeline run, if this is not given a datestamp is ued instead
#The Awk gets the N line from each text file in the input folder, i.e the individual stat files, the summary stats file is one line behind, 
#so this line number doesn't exist in the summary file
#The final sed removes the last line and creates the final csv file
# Author Steve Walsh June 2024
###############################################################################################################################################

# Function to handle incorrect arguments
function exit_with_bad_args {
    echo "Usage: bash lane_merger.bash optional args: --input <input path> --id <Run ID>"
    echo "Invalid arguments provided" >&2
    exit # this stops the terminal closing when run as source
}

INPUT=~/out
RUNID=RUNID="Stats_summary-$(date '+%Y-%m-%d-%R')"

#Set the possible input options
options=$(getopt -o '' -l input: -l id: -- "$@") || exit_with_bad_args

#Get the inputs
eval set -- "$options"
while true; do
    case "$1" in
        --input)
            shift
            INPUT="$1"
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

cd ${INPUT}/stats
echo \#Run ID\: ${RUNID} > summary_stats.txt
echo -n Sample_name, >> summary_stats.txt
awk 'FNR==3{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Forward_fastq, >> summary_stats.txt
awk 'FNR==4{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Reversqe_fastq, >> summary_stats.txt
awk 'FNR==5{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Aligned_reads, >> summary_stats.txt
awk 'FNR==6{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Aligned_percentage, >> summary_stats.txt
awk 'FNR==7{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Q30_Aligned_reads, >> summary_stats.txt
awk 'FNR==8{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Q30_Aligned_percentage, >> summary_stats.txt
awk 'FNR==9{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Q10_Aligned_reads, >> summary_stats.txt
awk 'FNR==10{printf $0 >> "summary_stats.txt"}' *
echo "" >> summary_stats.txt
echo -n Q10_Aligned_percentage, >> summary_stats.txt
awk 'FNR==11{printf $0 >> "summary_stats.txt"}' *
sed 's/.$//' summary_stats.txt >> ${RUNID}_summary_stats.csv
rm summary_stats.txt
