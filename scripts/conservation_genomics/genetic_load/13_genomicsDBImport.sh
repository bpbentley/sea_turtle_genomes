#!/bin/bash
#BSUB -J "GenomicsDBImport_CheMyd_exons[1-28]"
#BSUB -e ./logs/GenomicsDBImport/GenomicsDBImport_CheMyd_exons_%I.err
#BSUB -o ./logs/GenomicsDBImport/GenomicsDBImport_CheMyd_exons_%I.log
#BSUB -W 36:00
#BSUB -n 8
#BSUB -R rusage[mem=6000]
#BSUB -R span[hosts=1]
#BSUB -q long


module load gatk/4.1.8.1
module load htslib/1.9
module load samtools/1.4.1
module load python3/3.5.0

source /share/pkg/condas/2018-05-11/bin/activate && conda activate gatk_4.1.8.1

INDIR=/project/uma_lisa_komoroske/Blair/WGR_for_genome/snpEff/CheMyd/exons/
gatk --java-options "-Xmx40g -Xms40g" GenomicsDBImport \
      -V ${INDIR}/CheMyd_draft.g.vcf.gz \
      -V ${INDIR}/rCheMyd1.g.vcf.gz \
      -V ${INDIR}/SRR12153475.g.vcf.gz \
      -V ${INDIR}/SRR12153482.g.vcf.gz \
      -V ${INDIR}/SRR12153484.g.vcf.gz \
      -V ${INDIR}/SRR12153486.g.vcf.gz \
      --genomicsdb-workspace-path /project/uma_lisa_komoroske/Blair/WGR_for_genome/snpEff/CheMyd/exons/combined/SUPER_${LSB_JOBINDEX} \
      -L SUPER_${LSB_JOBINDEX}
      
conda deactivate && conda deactivate