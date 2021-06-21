#!/bin/bash
#BSUB -J rDerCor1_add_RG
#BSUB -e ./logs/rDerCor1_add_RG.err
#BSUB -o ./logs/rDerCor1_add_RG.log
#BSUB -W 24:00
#BSUB -n 16
#BSUB -R rusage[mem=4000]
#BSUB -R span[hosts=1]
#BSUB -q long

module load samtools/1.4.1
module load picard/2.23.3

BAMDIR=/project/uma_lisa_komoroske/Blair/rDerCor1/bams/BWA_rDerCor1_20210524
BAMIN=rDerCor1_20210524_SORT_DR.bam
BAMOUT=rDerCor1_20210524_SORT_DR_RG.bam

java -jar -Xms50G -Xmx60G /share/pkg/picard/2.23.3/picard.jar AddOrReplaceReadGroups \
 I=$BAMDIR/$BAMIN O=$BAMDIR/$BAMOUT \
 SORT_ORDER=coordinate RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=rDerCor1_20210524 CREATE_INDEX=True
 
samtools index $BAMDIR/$BAMOUT