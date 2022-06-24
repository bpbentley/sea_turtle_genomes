#!/bin/bash
#BSUB -J "combine_GVCFs_CheMyd[6]"
#BSUB -e ./logs/combine_GVCFs_CheMyd_%I.err
#BSUB -o ./logs/combine_GVCFs_CheMyd_%I.log
#BSUB -W 10:00
#BSUB -n 16
#BSUB -R rusage[mem=6000]
#BSUB -R span[hosts=1]
#BSUB -q long

DERCOR=/project/uma_lisa_komoroske/Blair/WGR_for_genome/GATK_out/CheMyd
SAMPLE=$(sed -n ${LSB_JOBINDEX}p ./WGR_greens.txt | cut -f1)
#SAMPLE=rCheMyd_genome
OUTDIR=./snpEff/CheMyd/exons

zcat ${DERCOR}/${SAMPLE}/exons/temp/${SAMPLE}_exons_01.g.vcf.gz | grep '#' > $OUTDIR/${SAMPLE}.g.vcf
for q in {01..28}; do
  zcat $DERCOR/${SAMPLE}/exons/temp/${SAMPLE}_exons_${q}.g.vcf.gz | grep -v '#' >> $OUTDIR/${SAMPLE}.g.vcf
done