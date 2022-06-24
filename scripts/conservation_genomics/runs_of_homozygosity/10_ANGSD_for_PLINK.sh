#!/bin/bash
#BSUB -J 'WGR_generate_PLINK_CheMyd[1-28]'
#BSUB -o ./logs/ANGSD_PLINK/WGR_CheMyd_ANGSD_PLINK.%I.log
#BSUB -e ./logs/ANGSD_PLINK/WGR_CheMyd_ANGSD_PLINK.%I.err
#BSUB -n 20
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=6000]
#BSUB -W 100:00
#BSUB -q long

module load samtools/1.4.1
module load angsd/0.933 
module load bedtools/2.29.2

REFDIR=/project/uma_lisa_komoroske/Blair/refs/rCheMyd1_20210528
REF=rCheMyd1.pri.cur.20210528.fasta
NUM=$(printf ${LSB_JOBINDEX})
CHR=SUPER_${NUM}

angsd -bam ./all_CheMyd_bams.txt -out ANGSD/CheMyd/CheMyd_${CHR} -doPlink 2 -doGeno -4 -doPost 1 -doMajorMinor 1 -GL 1 -doCounts 1 -doMaf 2 -postCutoff 0.95 -SNP_pval 1e-6 -uniqueonly 1 -remove_bads 1 -C 50 -baq 1 -r $CHR -ref $REFDIR/$REF -nThreads 20