#!/bin/bash

###############################################################################################################################################
############################## bulk rna-seq bash pipeline #####################################################################################
# This code will carryout a basic bulk rna-seq analysis pipeline for a pair of fastq files or pairs of fastq files fed from the sbatch script
# The pipeline can be initiated on the HPC using srun or locally with bash
# The pipeline assumes fastq file have been merged using the bash script on github
# The plan going forward is to implement Ahringer pipelines in Nextflow so this will not be developed beyond a basic workflow.
# Options include:
#      fastqid = The fastq file id, i.e. the start of the standard file namme without _merged_L1/2_001.fastq.gz
#      threads = Will multi-thread any process to this number
#      input = Change the path of the input fastq files, default is ~/data
#      id = Change the name of the output folder, the default is a datestamp
#      mergeID = If the file names have been merged differently the input can be changed here 'fastqid_<Add the flag here>_R1/R2_001.fastq.gz'
#      star_index = The location of the STAR index
#      kallisto_index = The location of the Kallisto index
# The scripts lacks logs and error handling
# Author Steve Walsh May 2024
################################################################################################################################################

#Set the defaults
outdir=~/out
fastq_dir=~/data/
BOWTIE_INDEX=~/references/built_genomes/star/c.elegans.latest
CHROM_SIZES=/mnt/home3/ahringer/index_files/genomes/c_elegans.PRJNA13758.WS285.genomic.chrom.sizes
THREADS=1
RUNID="PipelineRun-$(date '+%Y-%m-%d-%R')"
MERGEID=merged
base=null

# Function to handle incorrect arguments
function exit_with_bad_args {
    echo "Usage: bash lane_merger.bash optional args: --fastqid <fastq suffix> --sample_id <sample_id> --threads <number of threads> --input <input path> --id <Run ID>  --mergeID <merge ID> --bowtie_index"
    echo "Invalid arguments provided" >&2
    exit # this stops the terminal closing when run as source
}

#Set the possible input options
options=$(getopt -o '' -l fastqid: -l sample_id: -l threads: -l input: -l id: -l mergeID: -l bowtie_index: -- "$@") || exit_with_bad_args

#Get the inputs
eval set -- "$options"
while true; do
    case "$1" in
        --fastqid)
            shift
            FASTQ_ID="$1"
            ;;
        --sample_id)
            shift
            base="$1"
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

if [ ${base} = 'null' ]; then
    echo "No sample ID entered using fastq file name as ID"
    base=${FASTQ_ID}
else
    base=${FASTQ_ID}_${base}
fi

echo "${FASTQ_ID}${MERGEID}_R1_001.fastq.gz"
echo "${FASTQ_ID}${MERGEID}_R2_001.fastq.gz"

analysis_out_dir=${outdir}/${RUNID}
STATSFILE=${analysis_out_dir}/stats.csv

#Make the required output directories
mkdir ${analysis_out_dir}/${base}
mkdir ${analysis_out_dir}/${base}/fastq
mkdir ${analysis_out_dir}/${base}/trim_galore
trimmedfastq_dir=${analysis_out_dir}/${base}/trim_galore
mkdir ${analysis_out_dir}/${base}/star
mkdir ${analysis_out_dir}/${base}/fastq_screen
mkdir ${analysis_out_dir}/${base}/kallisto
cd ${analysis_out_dir}/${base}/fastq
cp $fastq_dir/${FASTQ_ID}${MERGEID}_R*_001.fastq.gz .

#Set up stats file
STATSFILE=${analysis_out_dir}/stats/stats-${base}.csv
echo \#Run ID,${RUNID} >> $STATSFILE
echo \# >> $STATSFILE
echo ${base}, >> $STATSFILE

#Gather fastq read numbers and add to stats file
R1count=$(( $(gunzip -c ${analysis_out_dir}/${base}/fastq/*R1_*.fastq.gz|wc -l)/4|bc ))
R2count=$(( $(gunzip -c ${analysis_out_dir}/${base}/fastq/*R2_*.fastq.gz|wc -l)/4|bc ))
echo ${R1count}, >> $STATSFILE
echo ${R2count}, >> $STATSFILE

#Carry out trimgalore (includes fastqc)
trim_galore --fastqc ${analysis_out_dir}/${base}/fastq/${FASTQ_ID}${MERGEID}_R1_001.fastq.gz ${analysis_out_dir}/${base}/fastq/${FASTQ_ID}${MERGEID}_R2_001.fastq.gz \
-o ${analysis_out_dir}/${base}/trim_galore \
-j ${THREADS}

#Carry out fastq screen
fastq_screen ${trimmedfastq_dir}/*.fq.gz  \
--outdir ${analysis_out_dir}/${base}/fastq_screen \
--threads ${THREADS}

#Carry out Bowtie alignment

echo "Carrying out Bowtie alignment"


#Filter to q30 reads
#samtools view -q 30 -b -h ${analysis_out_dir}/${base}/star/${base}_Aligned.sortedByCoord.out.bam > ${analysis_out_dir}/${base}/star/${base}_Aligned.sortedByCoord.out.q30.bam

#Index the bam file
samtools index  ${analysis_out_dir}/${base}/star/${base}_Aligned.sortedByCoord.out.bam
#samtools index  ${analysis_out_dir}/${base}/star/${base}_Aligned.sortedByCoord.out.q30.bam

#Add alignment stats to stats file
ALIGNEDREADS=$(samtools flagstat ${analysis_out_dir}/${base}/star/${base}_Aligned.sortedByCoord.out.bam)
Q30ALIGNEDREADS=$(samtools view -q 30 ${analysis_out_dir}/${base}/star/${base}_Aligned.sortedByCoord.out.bam | wc -l )
Q10ALIGNEDREADS=$(samtools view -q 10 ${analysis_out_dir}/${base}/star/${base}_Aligned.sortedByCoord.out.bam | wc -l )

ALIGNEDLIST=$(awk '{print $1;}' <<< "$ALIGNEDREADS")
Q30ALIGNEDREADSLIST=$(awk '{print $1;}' <<< "$Q30ALIGNEDREADS")
Q10ALIGNEDREADSLIST=$(awk '{print $1;}' <<< "$Q10ALIGNEDREADS")

ALIGNEDNUMBER=$(head -n 1 <<< $ALIGNEDLIST)
Q30ALIGNEDNUMBER=$(head -n 1 <<< $Q30ALIGNEDREADSLIST)
Q10ALIGNEDNUMBER=$(head -n 1 <<< $Q10ALIGNEDREADSLIST)
echo ${ALIGNEDNUMBER}, >> $STATSFILE
echo ${Q30ALIGNEDNUMBER}, >> $STATSFILE
echo ${Q10ALIGNEDREADS}, >> $STATSFILE

#Carry out MACS peak calling
