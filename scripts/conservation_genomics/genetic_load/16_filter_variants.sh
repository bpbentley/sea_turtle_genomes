#!/bin/bash
#BSUB -J filter_VCF_for_snpEff
#BSUB -o ./logs/filter_VCF_for_snpEff.log
#BSUB -e ./logs/filter_VCF_for_snpEff.err
#BSUB -n 8
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=5000]
#BSUB -W 05:00
#BSUB -q interactive
#BSUB -Is singularity

module load vcftools/0.1.16

VCFTOOLSIMG=/share/pkg/vcftools/0.1.16/vcftools-0.1.16.sif

for CHR in {9..28}; do

input_vcf=/project/uma_lisa_komoroske/Blair/WGR_for_genome/snpEff/DerCor/inputs/all_DerCor_SUPER_${CHR}.vcf
fileprefix=echo $CHR

######################################
### Filter VCF for use with snpEff ###
######################################

### 
singularity exec $VCFTOOLSIMG vcftools --vcf $input_vcf --min-meanDP 5 --max-meanDP 200 --recode --mac 1 --out ./snpEff/DerCor/inputs/SUPER_${fileprefix}

done