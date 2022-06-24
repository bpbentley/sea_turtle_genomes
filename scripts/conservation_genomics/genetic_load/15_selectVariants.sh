#!/bin/bash
#BSUB -J "DerCor_selectVariants[28]"
#BSUB -e ./logs/selectVariants/DerCor_selectVariants.%I.%J.err
#BSUB -o ./logs/selectVariants/DerCor_selectVariants.%I.%J.log
#BSUB -W 90:00
#BSUB -n 16
#BSUB -R rusage[mem=6000]
#BSUB -R span[hosts=1]
#BSUB -q long

# Load required modules
module load gatk/4.1.8.1
module load htslib/1.9
module load samtools/1.4.1
module load python3/3.5.0
module load htslib/1.9

source /share/pkg/condas/2018-05-11/bin/activate && conda activate gatk_4.1.8.1

# Adjust NAME to whatever prefix you want your output files to have
VCF=all_DerCor_SUPER_${LSB_JOBINDEX}.vcf # unique name needed if running multiple pipelines on different data, as files will be written to same directory, and will over-write if not unique.
TVCF=DerCor_SV_SUPER_${LSB_JOBINDEX}.vcf
CHR=SUPER_${LSB_JOBINDEX}

# Directories
INDIR=/project/uma_lisa_komoroske/Blair/WGR_for_genome/snpEff/DerCor/inputs
OUTDIR=/project/uma_lisa_komoroske/Blair/WGR_for_genome/snpEff/DerCor/inputs/selectVariants
REFDIR=/project/uma_lisa_komoroske/Blair/refs/rDerCor1_20210524

# Files
REFERENCE=rDerCor1.pri.cur.20210524.fasta  # reference must be indexed with samtools faidx, and data dictionary (samtools dict) first


gatk \
    SelectVariants \
    --java-options "-Xms60G -Xmx90G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    -R ${REFDIR}/${REFERENCE} \
    --select-type-to-include SNP \
    --select-type-to-include INDEL \
    -L ${CHR} \
    -V ${INDIR}/${VCF} \
    -O ${OUTDIR}/${TVCF}
    
conda deactivate && conda deactivate