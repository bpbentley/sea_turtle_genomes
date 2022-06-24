#!/bin/bash
#BSUB -J "CheMyd_WGR_align[5]"
#BSUB -o ./logs/CheMyd_WGR_align_%I.log
#BSUB -e ./logs/CheMyd_WGR_align_%I.err
#BSUB -q long
#BSUB -W 100:00
#BSUB -R rusage[mem=4000]
#BSUB -R span[hosts=1]
#BSUB -n 16

module load bwa/0.7.17
module load samtools/1.9
module load python3/3.8.2
module load fastqc/0.11.5
module load trimmomatic/0.39
module load picard/2.23.3

REFDIR=/project/uma_lisa_komoroske/Blair/refs/rCheMyd1_20210528
REF=rCheMyd1.pri.cur.20210528.fasta
FQDIR=/project/uma_lisa_komoroske/Blair/WGR_for_genome/fastq
SAMPLE=$(printf $(sed -n ${LSB_JOBINDEX}'p' ./WGR_greens.txt))
OUTDIR=./bams
METDIR=./metrics/alignment_metrics
METDIR2=./metrics/depth_metrics
TRIMDIR=./trimmed

# Run FastQC on the raw reads
#fastqc $FQDIR/${SAMPLE}_1.fastq.gz --outdir ./FastQC/raw_reads
#fastqc $FQDIR/${SAMPLE}_2.fastq.gz --outdir ./FastQC/raw_reads

# Run Trimmomatic to trim for quality
#java -jar /share/pkg/trimmomatic/0.39/trimmomatic-0.39.jar PE -phred33 -threads 16 $FQDIR/${SAMPLE}_1.fastq.gz $FQDIR/${SAMPLE}_2.fastq.gz \
#$TRIMDIR/paired/${SAMPLE}_1.paired.fastq.gz $TRIMDIR/unpaired/${SAMPLE}_1.unpaired.fastq.gz \
#$TRIMDIR/paired/${SAMPLE}_2.paired.fastq.gz $TRIMDIR/unpaired/${SAMPLE}_2.unpaired.fastq.gz \
#ILLUMINACLIP:/share/pkg/trimmomatic/0.39/adapters/NexteraPE-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:75

# Run FastQC on the trimmed reads
#fastqc $TRIMDIR/paired/${SAMPLE}_1.paired.fastq.gz --outdir ./FastQC/trimmed_reads
#fastqc $TRIMDIR/paired/${SAMPLE}_2.paired.fastq.gz --outdir ./FastQC/trimmed_reads

# Align trimmed reads to reference genome
bwa mem -t 16 $REFDIR/$REF $TRIMDIR/paired/${SAMPLE}_1.paired.fastq.gz $TRIMDIR/paired/${SAMPLE}_2.paired.fastq.gz > $OUTDIR/${SAMPLE}.sam

# Conver to BAM and delete SAM file
samtools view -S -b $OUTDIR/${SAMPLE}.sam > $OUTDIR/${SAMPLE}.bam
rm $OUTDIR/${SAMPLE}.sam

# Use SAMtools flagstat to get alignment statisctics:
samtools flagstat $OUTDIR/${SAMPLE}.bam > $METDIR/${SAMPLE}_alignment_stats.txt

# Sort and index file
samtools sort -o $OUTDIR/${SAMPLE}_SORT.bam $OUTDIR/${SAMPLE}.bam
samtools index $OUTDIR/${SAMPLE}_SORT.bam

# Check mean depth
java -jar -Xms50G -Xmx60G /share/pkg/picard/2.23.3/picard.jar CollectWgsMetrics \
       I=$OUTDIR/${SAMPLE}_SORT.bam \
       O=$METDIR2/$SAMPLE'_depth_metrics.txt' \
       R=$REFDIR/$REF