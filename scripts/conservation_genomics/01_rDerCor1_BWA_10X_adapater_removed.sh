#!/bin/bash
#BSUB -J rDerCor1_BWA_May2021_genome
#BSUB -o ./logs/rDerCor1_BWA_May2021_genome.%J.log
#BSUB -e ./logs/rDerCor1_BWA_May2021_genome.%J.err
#BSUB -q long
#BSUB -W 200:00
#BSUB -R rusage[mem=8000]
#BSUB -R span[hosts=1]
#BSUB -n 30

module load bwa/0.7.17
module load samtools/1.9
module load python3/3.8.2

REFDIR=/project/uma_lisa_komoroske/Blair/refs/rDerCor1_20210524
REF=rDerCor1.pri.cur.20210524.fasta.gz
FQDIR=../10X_removed_tags/
SAMPLE=rDerCor1
OUTDIR=./bams/BWA_rDerCor1_20210525
METDIR=/project/uma_lisa_komoroske/Blair/rDerCor1/metrics/alignment_metrics

# Concatenate the R1 and R2 reads:
#zcat $FQDIR/${SAMPLE}_S1_L001_R1_001.fastq.gz $FQDIR/${SAMPLE}_S1_L002_R1_001.fastq.gz > $FQDIR/${SAMPLE}_10XREM_S1_CAT_R1_001.fastq.gz
#zcat $FQDIR/${SAMPLE}_S1_L001_R2_001.fastq.gz $FQDIR/${SAMPLE}_S1_L002_R2_001.fastq.gz > $FQDIR/${SAMPLE}_10XREM_S1_CAT_R2_001.fastq.gz

# Index reference genome (one off)
bwa index $REFDIR/$REF

# Align reads to reference genome:
bwa mem -t 30 $REFDIR/$REF $FQDIR/${SAMPLE}_10XREM_S1_CAT_R1_001.fastq.gz $FQDIR/${SAMPLE}_10XREM_S1_CAT_R2_001.fastq.gz > $OUTDIR/${SAMPLE}_May2021_genome.sam

# Convert SAM to BAM file and sort:
samtools sort $OUTDIR/${SAMPLE}_May2021_genome.sam -o $OUTDIR/${SAMPLE}_May2021_genome_SORT.bam
samtools index $OUTDIR/${SAMPLE}_May2021_genome_SORT.bam

# Check alignment metrics:
samtools flagstat $OUTDIR/${SAMPLE}_May2021_genome_SORT.bam > $METDIR/${SAMPLE}_May2021_genome_alignment_stats.txt