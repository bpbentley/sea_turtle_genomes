#!/bin/bash
#BSUB -J rDerCor1_indel_realign
#BSUB -o ./logs/rDerCor1_indel_realign.log
#BSUB -e ./logs/rDerCor1_indel_realign.err
#BSUB -n 20
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=8000]
#BSUB -W 150:00
#BSUB -q long

module load gatk/3.5
module load htslib/1.9
module load samtools/1.4.1
module load python3/3.8.2
module load htslib/1.9
module load picard/2.23.3

BAMDIR=/project/uma_lisa_komoroske/Blair/rDerCor1/bams/BWA_rDerCor1_20210524
BAMIN=rDerCor1_20210524_SORT_DR_RG.bam
BAMOUT=rDerCor1_20210524_SORT_DR_RG_IR.bam
REFDIR=/project/uma_lisa_komoroske/Blair/refs/rDerCor1_20210524
REF=rDerCor1.pri.cur.20210524.fasta

#gunzip /project/uma_lisa_komoroske/DC_genome/rDerCor1_May2021_version/rDerCor1.pri.cur.20210524.fasta.gz
#ln -s /project/uma_lisa_komoroske/DC_genome/rDerCor1_May2021_version/rDerCor1.pri.cur.20210524.fasta $REFDIR

# First make sure the reference has been indexed
#samtools faidx $REFDIR/$REF

# Then create a sequence dictionary
java -jar -Xms120G -Xmx150G /share/pkg/picard/2.23.3/picard.jar CreateSequenceDictionary R=$REFDIR/$REF

## Realign around indels
## Create list of potential indels
## Note that RealignerTargetCreator was removed in GATK v >4.0.0
## Require version 3.5 for initial step
java -jar -Xms120G -Xmx150G /share/pkg/GATK/3.5/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R $REFDIR/$REF \
-I $BAMDIR/$BAMIN \
-o $BAMDIR'/rDerCor1_20210528.intervals' \
-drf BadMate


## Run the indel realigner tool
java -jar -Xms120G -Xmx150G /share/pkg/GATK/3.5/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R $REFDIR/$REF \
-I $BAMDIR/$BAMIN \
-targetIntervals $BAMDIR'/rDerCor1_20210528.intervals' \
--consensusDeterminationModel USE_READS  \
-o $BAMDIR/$BAMOUT