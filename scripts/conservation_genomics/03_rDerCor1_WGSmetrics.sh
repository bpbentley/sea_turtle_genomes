#!/bin/bash
#BSUB -J rDerCor1_WGS_metrics
#BSUB -o ./logs/rDerCor1_WGS_metrics.log
#BSUB -e ./logs/rDerCor1_WGS_metrics.err
#BSUB -q long
#BSUB -W 48:00
#BSUB -R rusage[mem=6000]
#BSUB -R span[hosts=1]
#BSUB -n 12

module load picard/2.23.3

REFDIR=/project/uma_lisa_komoroske/Blair/refs/rDerCor1_20210524
REF=rDerCor1.pri.cur.20210524.fasta.gz
BAMDIR=/project/uma_lisa_komoroske/Blair/rDerCor1/bams/BWA_rDerCor1_20210524
BAM=rDerCor1_20210524_SORT_DR.bam
METDIR=/project/uma_lisa_komoroske/Blair/rDerCor1/metrics/depth_metrics

java -Xms50G -Xmx60G  -jar /share/pkg/picard/2.23.3/picard.jar CollectWgsMetrics \
       I=$BAMDIR/$BAM \
       O=$METDIR/rDerCor1_20210524_depth_metrics.txt \
       R=$REFDIR/$REF
