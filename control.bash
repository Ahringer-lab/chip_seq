#!/bin/bash

###############################################################################################################################################
############################## bulk rna-seq bash pipeline #####################################################################################
# This code will carryout alignment of the control samples if required, the resulting bams can then be used across multiple treated samples
# This script can be used on its own if you just want to get bams for control samples
# This script would form part of a full pipeline
###############################################################################################################################################

#Set the defaults
outdir=~/out
fastq_dir=~/data/
BOWTIE_INDEX=/mnt/home3/ahringer/index_files/built_indexes/bwa
CHROM_SIZES=/mnt/home3/ahringer/index_files/genomes/c_elegans.PRJNA13758.WS285.genomic.chrom.sizes
THREADS=1
RUNID="PipelineRun-$(date '+%Y-%m-%d-%R')"
MERGEID=merged
WD="$(pwd)"

# Function to handle incorrect arguments
function exit_with_bad_args {
    echo "Usage: bash lane_merger.bash optional args: --fastqid <fastq suffix> --sample_id <sample_id> --threads <number of threads> --input <input path> --id <Run ID>  --mergeID <merge ID> --bowtieindex --callpeaks <true/false> --control <control for peak calling>  --controlbam <control bam file>" 
    echo "Invalid arguments provided" >&2
    exit # this stops the terminal closing when run as source
}

#Set the possible input options
options=$(getopt -o '' -l fastqid: -l sample_id: -l threads: -l input: -l id: -l mergeID: -l bowtieindex: -l "$@") || exit_with_bad_args

#Get the inputs
eval set -- "$options"
while true; do
    case "$1" in
        --fastqid)
            shift
            FASTQ_ID="$1"
            ;;
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
        --bowtieindex)
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

#Set up stats file
STATSFILECONTROL=${analysis_out_dir}/control/stats-${FASTQ_ID}.csv
echo \#Run ID,${RUNID} >> $STATSFILECONTROL
echo \# >> $STATSFILECONTROL
echo ${base}, >> $STATSFILECONTROL

#Gather fastq read numbers and add to stats file
R1count=$(( $(gunzip -c ${analysis_out_dir}/control/fastq/*R1_*.fastq.gz|wc -l)/4|bc ))
R2count=$(( $(gunzip -c ${analysis_out_dir}/control/fastq/*R2_*.fastq.gz|wc -l)/4|bc ))
echo ${R1count}, >> $STATSFILECONTROL
echo ${R2count}, >> $STATSFILECONTROL

#Carry out trimgalore (includes fastqc)
trim_galore --fastqc ${analysis_out_dir}/control/${FASTQ_ID}${MERGEID}_R1_001.fastq.gz ${analysis_out_dir}/control/fastq/${FASTQ_ID}${MERGEID}_R2_001.fastq.gz \
-o ${analysis_out_dir}/${base}/trim_galore \
-j ${THREADS}

#Carry out fastq screen
fastq_screen ${trimmedfastq_dir}/*.fq.gz  \
--outdir ${analysis_out_dir}/${base}/fastq_screen \
--threads ${THREADS}

#Carry out BWA alignment
echo "Carrying out BWA alignment"
bwa mem -t ${THREADS} ${trimmedfastq_dir}/${FASTQ_ID}${MERGEID}_R*_001_trimmed.fq.gz > ${analysis_out_dir}/control/${base}.sam

#Convert to bam, not currently downsampling to q10 by default, add as option?
samtools view -@ ${THREADS} -b -h ${analysis_out_dir}/${base}/control/${base}.sam > ${analysis_out_dir}/control/${base}.bam
#samtools view -@ ${THREADS} -q 10 -b -h ${analysis_out_dir}/${base}/bwa/${base}.sam > ${analysis_out_dir}/${base}/bwa/${base}.q10.bam

#Sort the bam file
samtools sort -@ ${THREADS} ${analysis_out_dir}/${base}/control/${base}.bam ${analysis_out_dir}/${base}/control/${base}.sorted.bam
#samtools sort -@ ${THREADS} ${analysis_out_dir}/${base}/bwa/${base}.q10.bam ${analysis_out_dir}/${base}/bwa/${base}.q10.sorted.bam

#Index the bam files
samtools index  ${analysis_out_dir}/${base}/control/${base}.sorted.bam
#samtools index  ${analysis_out_dir}/${base}/bwa/${base}.sorted.bam

#Clean up files
rm ${analysis_out_dir}/${base}/control/${base}.sam
rm ${analysis_out_dir}/${base}/control/${base}.bam
#rm ${analysis_out_dir}/${base}/bwa/${base}.q10.sam


#Add alignment stats to stats file
ALIGNEDREADS=$(samtools flagstat ${analysis_out_dir}/${base}/control/${base}.sorted.bam)
Q30ALIGNEDREADS=$(samtools view -q 30 ${analysis_out_dir}/${base}/control/${base}.sorted.bam | wc -l )
Q10ALIGNEDREADS=$(samtools view -q 10 ${analysis_out_dir}/${base}/control/${base}.sorted.bam | wc -l )

ALIGNEDLIST=$(awk '{print $1;}' <<< "$ALIGNEDREADS")
Q30ALIGNEDREADSLIST=$(awk '{print $1;}' <<< "$Q30ALIGNEDREADS")
Q10ALIGNEDREADSLIST=$(awk '{print $1;}' <<< "$Q10ALIGNEDREADS")

ALIGNEDNUMBER=$(head -n 1 <<< $ALIGNEDLIST)
Q30ALIGNEDNUMBER=$(head -n 1 <<< $Q30ALIGNEDREADSLIST)
Q10ALIGNEDNUMBER=$(head -n 1 <<< $Q10ALIGNEDREADSLIST)
echo ${ALIGNEDNUMBER}, >> $STATSFILECONTROL
echo ${Q30ALIGNEDNUMBER}, >> $STATSFILECONTROL
echo ${Q10ALIGNEDREADS}, >> $STATSFILECONTROL