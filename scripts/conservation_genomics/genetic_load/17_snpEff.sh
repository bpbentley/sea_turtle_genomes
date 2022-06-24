#!/bin/bash
#BSUB -J "CheMyd_snpEff_sep[1-5]"
#BSUB -o ./logs/snpEff/CheMyd_snpEff_sep.%I.log
#BSUB -e ./logs/snpEff/CheMyd_snpEff_sep.%I.err
#BSUB -n 20
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=6000]
#BSUB -W 24:00
#BSUB -q long

module load snpEff_snpSift/4.3T
module load gatk/4.1.8.1
module load htslib/1.9
module load samtools/1.4.1
module load python3/3.5.0
module load vcftools/0.1.16

DERCOR=/project/uma_lisa_komoroske/Blair/WGR_for_genome/GATK_out/CheMyd
OUTDIR=/project/uma_lisa_komoroske/Blair/WGR_for_genome/snpEff/CheMyd/exons/separate
snpEff=/share/pkg/snpEff_snpSift/4.3T/snpEff
config=/project/uma_lisa_komoroske/Blair/refs/snpEff/snpEff.config
DIR=/project/uma_lisa_komoroske/Blair/WGR_for_genome/snpEff/CheMyd/exons/separate
SAMPLE=$(sed -n ${LSB_JOBINDEX}p ./WGR_greens.txt | cut -f1)
REFDIR=/project/uma_lisa_komoroske/Blair/refs/rCheMyd1_20210528
REF=rCheMyd1.pri.cur.20210528.fasta
VCFTOOLSIMG=/share/pkg/vcftools/0.1.16/vcftools-0.1.16.sif


#zcat ${DERCOR}/${SAMPLE}/exons/${SAMPLE}_01.Filter.vcf.gz | grep '#' > $OUTDIR/${SAMPLE}.Filter.vcf
#for q in {01..28}; do
#  zcat $DERCOR/${SAMPLE}/exons/${SAMPLE}_${q}.Filter.vcf.gz | grep -v '#' >> $OUTDIR/${SAMPLE}.Filter.vcf
#done

#source /share/pkg/condas/2018-05-11/bin/activate && conda activate gatk_4.1.8.1

#gatk \
#    SelectVariants \
#    --java-options "-Xmx80G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
#    -R ${REFDIR}/${REF} \
#    --select-type-to-include SNP \
#    --select-type-to-include INDEL \
#    -V $OUTDIR/${SAMPLE}.Filter.vcf \
#    -O $OUTDIR/${SAMPLE}.Variants.vcf \
    
#conda deactivate && conda deactivate


#singularity exec $VCFTOOLSIMG vcftools --gzvcf $OUTDIR/${SAMPLE}.Variants.vcf --remove-filtered-all --recode --out $OUTDIR/${SAMPLE}.Variants.recode.vcf

java -Xmx80G -jar $snpEff/snpEff.jar eff rCheMyd1 -c ${config} -s $DIR/${SAMPLE}_exons \
$OUTDIR/${SAMPLE}.Variants.recode.vcf.recode.vcf > $DIR/${SAMPLE}.snpEff.vcf
