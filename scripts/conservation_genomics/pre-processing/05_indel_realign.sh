#!/bin/bash
#BSUB -J "04_WGR_CheMyd_indel_realign[1-4]"
#BSUB -o ./logs/indel_realign/greens/CheMyd_indel_realign_%I.log
#BSUB -e ./logs/indel_realign/greens/CheMyd_indel_realign_%I.err
#BSUB -q long
#BSUB -W 36:00
#BSUB -R rusage[mem=4000]
#BSUB -R span[hosts=1]
#BSUB -n 20

module load gatk/3.5
module load htslib/1.9
module load samtools/1.4.1
module load python3/3.8.2
module load htslib/1.9
module load picard/2.23.3

SAMPLE=$(sed -n ${LSB_JOBINDEX}'p' ./WGR_greens.txt)
REFDIR=/project/uma_lisa_komoroske/Blair/refs/rCheMyd1_20210528
REF=rCheMyd1.pri.cur.20210528.fasta
BAMDIR=/project/uma_lisa_komoroske/Blair/WGR_for_genome/bams

BAMIN=${SAMPLE}_SORT_DR_RG.bam
BAMOUT=${SAMPLE}_SORT_DR_RG_IR.bam

## Realign around indels
## Create list of potential indels
## Note that RealignerTargetCreator was removed in GATK v >4.0.0
## Require version 3.5 for initial step
java -jar -Xms120G -Xmx150G /share/pkg/GATK/3.5/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R $REFDIR/$REF \
-I $BAMDIR/$BAMIN \
-o $BAMDIR/rCheMyd1_20210528_${SAMPLE}.intervals \
-drf BadMate


## Run the indel realigner tool
java -jar -Xms120G -Xmx150G /share/pkg/GATK/3.5/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R $REFDIR/$REF \
-I $BAMDIR/$BAMIN \
-targetIntervals $BAMDIR/rCheMyd1_20210528_${SAMPLE}.intervals \
--consensusDeterminationModel USE_READS  \
-o $BAMDIR/$BAMOUT