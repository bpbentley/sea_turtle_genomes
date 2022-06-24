#!/bin/bash
#BSUB -J rCheMyd1_proc10x
#BSUB -o ./logs/rCheMyd1_proc10X.%I.log
#BSUB -e ./logs/rCheMyd1_proc10X.%I.err
#BSUB -q long
#BSUB -W 24:00
#BSUB -R rusage[mem=4000]
#BSUB -R span[hosts=1]
#BSUB -n 4

module load python/2.7.9

inputdir=./fastq
outdir=./fastq

./scripts/process_10xReads.py \
-1 ${inputdir}/rCheMyd1_R1.fastq.gz \
-2 ${inputdir}/rCheMyd1_R2.fastq.gz \
-o ${outdir}/rCheMyd1 -a