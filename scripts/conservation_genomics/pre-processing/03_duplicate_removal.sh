#!/bin/bash
#BSUB -J "WGR_CheMyd_duplicate_removal[6]"
#BSUB -o ./logs/dupe_remove/CheMyd_duplicate_removal_%I.log
#BSUB -e ./logs/dupe_remove/CheMyd_duplicate_removal_%I.err
#BSUB -q long
#BSUB -W 48:00
#BSUB -R rusage[mem=4000]
#BSUB -R span[hosts=1]
#BSUB -n 20

module load samtools/1.4.1
module load picard/2.23.3

BAMDIR=./bams
METDIR=./metrics/dupe_metrics
METDIR2=./metrics/depth_metrics/postdupe
SAMPLE=$(printf $(sed -n ${LSB_JOBINDEX}'p' ./WGR_greens.txt))
REFDIR=/project/uma_lisa_komoroske/Blair/refs/rCheMyd1_20210528
REF=rCheMyd1.pri.cur.20210528.fasta

java -Xms30G -Xmx60G -jar /share/pkg/picard/2.23.3/picard.jar MarkDuplicates \
      I=$BAMDIR/${SAMPLE}_SORT.bam \
      O=$BAMDIR/${SAMPLE}_SORT_DR.bam \
      M=$METDIR/${SAMPLE}_duplicate_metrics.txt \
      ASSUME_SORT_ORDER=coordinate \
      VALIDATION_STRINGENCY=SILENT \
      REMOVE_DUPLICATES=TRUE \
      MAX_RECORDS_IN_RAM=500000
      
samtools index $BAMDIR/${SAMPLE}_SORT_DR.bam

java -jar -Xms50G -Xmx60G /share/pkg/picard/2.23.3/picard.jar CollectWgsMetrics \
       I=$BAMDIR/${SAMPLE}_SORT_DR.bam \
       O=$METDIR2/$SAMPLE'_depth_metrics_post_DR.txt' \
       R=$REFDIR/$REF