#!/bin/bash
#BSUB -J "CheMyd_masked[1-28]"
#BSUB -e ./logs/GATK/masked/CheMyd_GATK_masked.%I.%J.err
#BSUB -o ./logs/GATK/masked/CheMyd_GATK_masked.%I.%J.log
#BSUB -W 90:00
#BSUB -n 20
#BSUB -R rusage[mem=6000]
#BSUB -R span[hosts=1]
#BSUB -q long

# Load required modules
module load gatk/4.1.8.1
module load htslib/1.9
module load samtools/1.4.1
module load python3/3.5.0
module load htslib/1.9

source /share/pkg/condas/2018-05-11/bin/activate && conda activate gatk_4.1.8.1

###############################################################################
# These values will need to be adjusted to fit your environment as well
# as on a per analysis run basis

# NOTE: Important!!!! Need to set the min and max depth for the genome ID in the python script filterVCF_20200723.py
	#   Good idea to try to modify this script to update the python script prior to running. 

# Adjust NAME to whatever prefix you want your output files to have
NAME=$(sed -n 5p ./WGR_greens.txt | cut -f1) # unique name needed if running multiple pipelines on different data, as files will be written to same directory, and will over-write if not unique.

# Directories
SCRIPTDIR=/project/uma_lisa_komoroske/Blair/rCheMyd1/scripts/gatk_in
BAMDIR=/project/uma_lisa_komoroske/Blair/WGR_for_genome/bams
OUTDIR=/project/uma_lisa_komoroske/Blair/WGR_for_genome/GATK_out/CheMyd/${NAME}/masked
TEMPDIR=/project/uma_lisa_komoroske/Blair/WGR_for_genome/GATK_out/CheMyd/${NAME}/masked/temp
REFDIR=/project/uma_lisa_komoroske/Blair/refs/rCheMyd1_20210528

# Files
REFERENCE=rCheMyd1.pri.cur.20210528.fasta  # reference must be indexed with samtools faidx, and data dictionary (samtools dict) first
SCAFFOLDLIST=rCheMyd1.pri.cur.20210528.scaffold.list.txt # list must have ">" removed from before each scaffold! SEE README FILE FOR HOW TO GENERATE THE SCAFFOLDLIST FROM THE FASTA FILE
CHRLENGTHS=${REFDIR}/rCheMyd1.pri.cur.20210528.scaffold.lengths.txt # SEE README FILE FOR HOW TO GENERATE THE CHROMOSOME LENGTHS FILE FROM THE FASTA INDEX FILE
BAM=${NAME}_SORT_DR_RG.bam # samtools index bam file first!

# scripts
PSCRIPT=${SCRIPTDIR}/filterVCF_010920.py # python script for filtering vcf file
WSCRIPT=${SCRIPTDIR}/slidingWindowHet_010920.py # python script for creating sliding windows and counting heterozygotes per window. Specify window size and step size:

# values
WINSIZE=100000
STEPSIZE=100000
# minimum and maximum depth of coverage in BAM file (usually set to 1/3x and 2x average depth of coverage; assumes one sample per BAM file). 
MINDEPTH=$(sed -n 5p ./WGR_greens.txt | cut -f3)
MAXDEPTH=$(sed -n 5p ./WGR_greens.txt | cut -f4)

###############################################################################
# Get working scaffold based on array number
NUM=$(printf %02d ${LSB_JOBINDEX})
CHR=$(head -n ${NUM} ${REFDIR}/${SCAFFOLDLIST} | tail -n 1)
# NC_045784.1 Phocoena sinus isolate mPhoSin1 chromosome X, mPhoSin1.pri, whole genome shotgun sequence
CHR=$(echo ${CHR} | awk -F " " '{ print $1 }')
echo ${CHR}

# Output Files
MYLOG=${OUTDIR}/${NAME}_gatk.${NUM}.log
GVCF=${TEMPDIR}/${NAME}_${NUM}.g.vcf.gz
VCF=${TEMPDIR}/${NAME}_${NUM}.vcf.gz
TVCF=${TEMPDIR}/${NAME}_${NUM}.TrimAlt.vcf.gz
AVCF=${TEMPDIR}/${NAME}_${NUM}.AddAnnot.vcf.gz
FVCF=${OUTDIR}/${NAME}_${NUM}.Filter.vcf.gz

###############################################################################

# Start pipeline
# HaplotypeCaller
echo -e "[$(date "+%Y-%m-%d %T")] Starting HaplotypeCaller" >> ${MYLOG}
gatk \
    HaplotypeCaller \
    --java-options "-Xms60G -Xmx100G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    -R ${REFDIR}/${REFERENCE} \
    -ERC BP_RESOLUTION \
    -mbq 20 \
    -L $REFDIR/masked_windows/${CHR}_masked.bed \
    -I ${BAMDIR}/${BAM} \
    -O ${GVCF} \
    --output-mode EMIT_ALL_ACTIVE_SITES \
    >> ${MYLOG} 2>&1

#############
# GenotypeGVCFs
# GATK3->GATK4 Update
# --include-non-variant-sites replaces -allSites
# --standard-min-confidence-threshold-for-calling replaces -stand_call_conf
echo -e "[$(date "+%Y-%m-%d %T")] Starting GenotypeGVCFs" >> ${MYLOG}
gatk \
    GenotypeGVCFs \
    --java-options "-Xms60G -Xmx100G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    -R ${REFDIR}/${REFERENCE} \
    --include-non-variant-sites \
    --standard-min-confidence-threshold-for-calling 0 \
    -L $REFDIR/masked_windows/${CHR}_masked.bed \
    -V ${GVCF} \
    -O ${VCF} \
    >> ${MYLOG} 2>&1

#############
# SelectVariants
# --remove-unused-alternates replaces -trimAlternates
echo -e "[$(date "+%Y-%m-%d %T")] Starting SelectVariants" >> ${MYLOG}
gatk \
    SelectVariants \
    --java-options "-Xms60G -Xmx100G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    -R ${REFDIR}/${REFERENCE} \
    --remove-unused-alternates \
    -L $REFDIR/masked_windows/${CHR}_masked.bed \
    -V ${VCF} \
    -O ${TVCF} \
    >> ${MYLOG} 2>&1

#############
# VariantAnnotator
# VariantType removed from GATK4 per
# https://gatkforums.broadinstitute.org/gatk/discussion/13500/a-varianttype-annotation-still-available-in-gatk4-haplotypecaller
# This is not functioning properly, and ends up changing the coverage per allele to 0,0 for all heterozygotes, so that there are none that pass filter in the next step.  Skip for now and run the final filter on the TrimAlt.vcf.gz files

# echo -e "[$(date "+%Y-%m-%d %T")] Starting VariantAnnotator" >> ${MYLOG}
# gatk \
#     VariantAnnotator \
#     --java-options "-Xms60G -Xmx100G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
#     -R ${REFDIR}/${REFERENCE} \
#     -G StandardAnnotation \
#     -L ${CHR} \
#     -V ${TVCF} \
#     -O ${AVCF} \
#     >> ${MYLOG} 2>&1

# python script to filter VCF (from AddAnnot.vcf.gz files)
# echo -e "[$(date "+%Y-%m-%d %T")] Starting Python Script" >> ${MYLOG}
# python3 ${PSCRIPT} ${AVCF} | bgzip > ${FVCF}
# tabix -p vcf ${FVCF}

#############
# Filter VCF
#python script to filter VCF (from TrimAlt.vcf.gz files)
echo -e "[$(date "+%Y-%m-%d %T")] Starting Python Script" >> ${MYLOG}
python3 ${PSCRIPT} ${TVCF} ${MINDEPTH} ${MAXDEPTH} | bgzip > ${FVCF}
tabix -p vcf ${FVCF}


echo -e "[$(date "+%Y-%m-%d %T")] Finished Pipeline" >> ${MYLOG}
wait

##############################################################################
# WinHet script
# uses python, and module "pysam" for reading, manipulating and writing genomic data sets.
# to install pysam, source your python venv
source /home/bb89a/python3/bin/activate
# pip install pysam

CHR=$(head -n ${NUM} ${CHRLENGTHS} | tail -n 1)

python ${WSCRIPT} ${FVCF} ${WINSIZE} ${STEPSIZE} ${CHR} ${CHRLENGTHS}

deactivate
conda deactivate && conda deactivate