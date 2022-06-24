#!/bin/bash
#BSUB -J "03_WGR_CheMyd_add_RG[1-5]"
#BSUB -o ./logs/add_RG/greens/CheMyd_add_RG_%I.log
#BSUB -e ./logs/add_RG/greens/_add_RG_%I.err
#BSUB -q long
#BSUB -W 36:00
#BSUB -R rusage[mem=4000]
#BSUB -R span[hosts=1]
#BSUB -n 20

module load samtools/1.4.1
module load picard/2.23.3

BAMDIR=/project/uma_lisa_komoroske/Blair/WGR_for_genome/bams
SAMPLE=$(sed -n ${LSB_JOBINDEX}'p' ./WGR_greens.txt)
REFDIR=/project/uma_lisa_komoroske/Blair/refs/rCheMyd1_20210528
REF=rCheMyd1.pri.cur.20210528.fasta

BAMIN=${SAMPLE}_SORT_DR.bam
BAMOUT=${SAMPLE}_SORT_DR_RG.bam

java -jar -Xms50G -Xmx60G /share/pkg/picard/2.23.3/picard.jar AddOrReplaceReadGroups \
 I=$BAMDIR/$BAMIN O=$BAMDIR/$BAMOUT \
 SORT_ORDER=coordinate RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${SAMPLE} CREATE_INDEX=True
 
samtools index $BAMDIR/$BAMOUT