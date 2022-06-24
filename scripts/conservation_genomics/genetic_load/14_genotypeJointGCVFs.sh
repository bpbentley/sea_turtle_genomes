#!/bin/bash
#BSUB -J "Genotype_jointVCF_CheMyd_exons[1-28]"
#BSUB -e ./logs/jointGVCF/jointGVCF_CheMyd_exons.%I.err
#BSUB -o ./logs/jointGVCF/jointGVCF_CheMyd_exons.%I.log
#BSUB -W 72:00
#BSUB -n 8
#BSUB -R rusage[mem=6000]
#BSUB -R span[hosts=1]
#BSUB -q long


module load gatk/4.1.8.1
module load htslib/1.9
module load samtools/1.4.1
module load python3/3.5.0

REFDIR=/project/uma_lisa_komoroske/Blair/refs/rCheMyd1_20210528
REF=rCheMyd1.pri.cur.20210528.fasta
OUTDIR=/project/uma_lisa_komoroske/Blair/WGR_for_genome/snpEff/CheMyd/exons/vcf

source /share/pkg/condas/2018-05-11/bin/activate && conda activate gatk_4.1.8.1

gatk --java-options "-Xmx40g -Xms40g" GenotypeGVCFs \
      -R ${REFDIR}/${REF} \
      -V gendb:///project/uma_lisa_komoroske/Blair/WGR_for_genome/snpEff/CheMyd/exons/combined/SUPER_${LSB_JOBINDEX} \
      -O ${OUTDIR}/CheMyd_exons_SUPER_${LSB_JOBINDEX}.vcf
      
conda deactivate && conda deactivate