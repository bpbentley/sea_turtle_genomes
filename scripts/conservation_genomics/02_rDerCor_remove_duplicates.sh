#!/bin/bash
#BSUB -J rDerCor1_20210524_dupe_remove
#BSUB -o ./logs/rDerCor1_20210524_dupe_remove.log
#BSUB -e ./logs/rDerCor1_20210524_dupe_remove.err
#BSUB -q long
#BSUB -W 48:00
#BSUB -R rusage[mem=5000]
#BSUB -R span[hosts=1]
#BSUB -n 20

module load samtools/1.4.1
module load picard/2.23.3

REFDIR=/project/uma_lisa_komoroske/Blair/refs/rDerCor1_20210524
REF=rDerCor1.pri.cur.20210524.fasta.gz
BAMDIR=/project/uma_lisa_komoroske/Blair/rDerCor1/bams/BWA_rDerCor1_20210524
BAM=rDerCor1_May2021_genome_SORT.bam
OUTBAM=rDerCor1_20210524_SORT_DR.bam
METDIR=/project/uma_lisa_komoroske/Blair/rDerCor1/metrics/dupe_metrics

# Ensure BAM file is sorted and indexed.
# Run Picard Tools to identify and remove duplicate reads:
#java -Xms60G -Xmx80G -jar /share/pkg/picard/2.23.3/picard.jar MarkDuplicates \
#      I=$BAMDIR/$BAM \
#      O=$BAMDIR/$OUTBAM \
#      M=$METDIR/rDerCor1_20210525_dupe_metrics.txt \
#      ASSUME_SORT_ORDER=coordinate \
#      REMOVE_DUPLICATES=TRUE MAX_RECORDS_IN_RAM=350000
      
samtools index $BAMDIR/$OUTBAM