#!/bin/bash
#BSUB -J "CheMyd_PSMC[1-5]"
#BSUB -o ./logs/PSMC/CheMyd_PSMC.%I.log
#BSUB -e ./logs/PSMC/CheMyd_PSMC.%I.err
#BSUB -q long
#BSUB -W 200:00
#BSUB -R rusage[mem=8000]
#BSUB -R span[hosts=1]
#BSUB -n 30

module load htslib/1.9
module load samtools/1.9
module load python3/3.5.0
module load htslib/1.9
module load anaconda3/2019.03
module load bedtools/2.29.2
module load gcc/8.1.0
module load bcftools/1.9
module load R/3.6.1
module load tabix/0.2.6
module load psmc/0.6.5
module load gnuplot/5.2.0

## Scripts required from Harvi:
# mpileup.vcf.v2.self.arr
# allelebalance.mpileupvcfs.arr
# allelebalance.mpileupvcfs.R
# filterout.ab.arr
# filterout_maskedpos.arr
# psmc_plot.pl


# 1) preprocessing the bams for input into psmc - ie call & filter SNPs

# 1. filter out the short fragments - ultrashort alignments from bwa mem
# 2. generate vcf
	#	with bcftools mpileup require mapping & base quality of 30  - # samtools view, bcftools mpileup, bcftools call

# amend script:
#/scratch/devel/hpawar/turtles/data/proc10xG/scripts/psmc/troubleshoot/mpileup.vcf.v2.self.arr

# new genome references
#rCheMyd1.pri.cur.20210528.fasta
#-----------------------------------------------------------------------------------------------------------------------
#Fri 13 Aug 2021 16:41:11 CEST - for bams mapped to final assembly
REFDIR=/project/uma_lisa_komoroske/Blair/refs/rCheMyd1_20210528
OUTDIR=/project/uma_lisa_komoroske/Blair/WGR_for_genome/PSMC/CheMyd
BAMDIR=/project/uma_lisa_komoroske/Blair/WGR_for_genome/bams
n=$(sed -n ${LSB_JOBINDEX}p ./WGR_greens.txt | cut -f1)
BAM=${n}_SORT_DR_RG_IR.bam
REF=rCheMyd1.pri.cur.20210528.fasta

samtools view -h $BAMDIR/$BAM | awk 'substr($0,1,1)=="@" || ($9>= 50 && $9<=5000) || ($9<=-50 && $9>=-5000)' | bcftools mpileup -Q 30 -q 30  -Ou -a "FORMAT/AD,FORMAT/DP,INFO/AD" -f $REFDIR/$REF - | bcftools call -c - > $OUTDIR/"$n".vcf

#-----------------------------------------------------------------------------------------------------------------------

# 3. assess allele balance of the vcfs generated above 
# using 
#/scratch/devel/hpawar/turtles/data/proc10xG/scripts/psmc/troubleshoot/allelebalance.mpileupvcfs.arr # calls:
#/scratch/devel/hpawar/turtles/data/proc10xG/scripts/psmc/troubleshoot/allelebalance.mpileupvcfs.R
# ie extract positions with AB<0.25, > 0.75
m=$(echo $n)

bgzip $OUTDIR/"$m".vcf > $OUTDIR/"$m".filt.vcf.gz
tabix $OUTDIR/"$m".filt.vcf.gz

mv $OUTDIR/"$m".vcf.gz $OUTDIR/"$m".filt.vcf.gz 
gunzip $OUTDIR/"$m".filt.vcf.gz

Rscript --vanilla /project/uma_lisa_komoroske/Blair/rCheMyd1/scripts/allelebalance.mpileupvcfs.R $m $OUTDIR

# formatting - AB<0.25, > 0.75
cat $OUTDIR/${m}_ab_pos.bed | tr -s ' ' | cut -d ' ' -f2,3,4 | sed 's/"//g' | sed '1d' | awk '{print $1"\t"$2"\t"$3}' > $OUTDIR/${m}_ab_pos2.bed
#-----------------------------------------------------------------------------------------------------------------------

# 4. filter out pos with AB<0.25, > 0.75 - /scratch/devel/hpawar/turtles/data/proc10xG/scripts/psmc/troubleshoot/filterout.ab.arr

bedtools subtract -a $OUTDIR/"$m".filt.vcf -b $OUTDIR/"$m"_ab_pos2.bed > $OUTDIR/"$m"_filt_ab.vcf.gz
#-----------------------------------------------------------------------------------------------------------------------

# formatting - to add the header back into _ab vcf 
#rCheMyd1_self_filt_ab.vcf.gz
touch $OUTDIR/header_${n}.hr

bcftools view -h  $OUTDIR/${m}.filt.vcf > $OUTDIR/header_${m}.hr

cat $OUTDIR/header_${m}.hr $OUTDIR/${m}_filt_ab.vcf.gz > $OUTDIR/${m}_filt_ab1.vcf.gz

#-----------------------------------------------------------------------------------------------------------------------

# 4.1 - filter by the masked positions - /scratch/devel/hpawar/turtles/data/proc10xG/scripts/psmc/troubleshoot/filterout_maskedpos.arr
# files contain the positions to be retained
#/scratch/devel/hpawar/turtles/data/bb_bams/finalassembly/rCheMyd1.pri.cur.20210528.masked.bed

#inputdir="/scratch/devel/hpawar/turtles/data/bb_bams/finalassembly/205.172.168.25/uma_pLKg561a/vcf_mpileup"
#maskedpos=$REFDIR/"rCheMyd1.pri.cur.20210528.masked.bed"

#bedtools intersect -header -a $OUTDIR/"$m"_filt_ab1.vcf.gz -b "$maskedpos"  > $OUTDIR/"$m"_filt_ab2.vcf.gz

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

# 5. vcfutils step - /scratch/devel/hpawar/turtles/data/proc10xG/scripts/psmc/troubleshoot/vcfutils.arr
# for bams mapped to final assembly version & additional filtering by masked pos
#inputdir="/scratch/devel/hpawar/turtles/data/bb_bams/finalassembly/205.172.168.25/uma_pLKg561a/vcf_mpileup"
# rCheMyd1_final_filt_ab2.vcf.gz
# rDerCor1_final_filt_ab2.vcf.gz

#if [[ "$ID" -eq 1 ]]; then m="rCheMyd1_final"; elif [[ "$ID" -eq 2 ]]; then m="rDerCor1_final"; else m="x";  fi;

# for green - # take 18-109
# for leatherback - take 30-181

minc=$(sed -n ${LSB_JOBINDEX}p ./WGR_greens.txt | cut -f3)
maxc=$(sed -n ${LSB_JOBINDEX}p ./WGR_greens.txt | cut -f4)

bcftools view -t SUPER_1,SUPER_2,SUPER_3,SUPER_4,SUPER_5,SUPER_6,SUPER_7,SUPER_8,SUPER_9,SUPER_10 $OUTDIR/"$m"_filt_ab1.vcf.gz -Oz -o $OUTDIR/"$m".10.vcf.gz 
tabix $OUTDIR/"$m".10.vcf.gz 

#-----------------------------------------------------------------------------------------------------------------------
#outdir="/scratch/devel/hpawar/turtles/psmc/proc10xG/extrafilt/pre_psmc"
#inputdir="/scratch/devel/hpawar/turtles/data/bb_bams/finalassembly/205.172.168.25/uma_pLKg561a/vcf_mpileup"
#if [[ "$ID" -eq 1 ]]; then m="rCheMyd1_final"; elif [[ "$ID" -eq 2 ]]; then m="rDerCor1_final"; else m="x";  fi;
# for green - # take 18-109
# for leatherback - take 30-181
#if [[ "$ID" -eq 1 ]]; then minc=18; elif [[ "$ID" -eq 2 ]]; then minc=30;  else minc=0;  fi;
#if [[ "$ID" -eq 1 ]]; then maxc=109; elif [[ "$ID" -eq 2 ]]; then maxc=181;  else maxc=0;  fi;
#module load BCFTOOLS/1.12
bcftools view $OUTDIR/"$m".10.vcf.gz  | /share/pkg/bcftools/1.9/bin/vcfutils.pl vcf2fq -d "$minc" -D "$maxc" -Q 30  | gzip > $OUTDIR/"$m".10supers."$minc"."$maxc".fq
#-----------------------------------------------------------------------------------------------------------------------

# 6. Run PSMC - /scratch/devel/hpawar/turtles/data/proc10xG/scripts/psmc/troubleshoot/psmc.step2.25aug.self.arr

# convert input fastq to input format required for PSMC (.psmcfa)
/share/pkg/psmc/0.6.5/utils/fq2psmcfa $OUTDIR/"$m".10supers."$minc"."$maxc".fq > $OUTDIR/"$m".10supers."$minc"."$maxc".psmcfa

# split for bootstrapping
/share/pkg/psmc/0.6.5/utils/splitfa $OUTDIR/"$m".10supers."$minc"."$maxc".psmcfa > $OUTDIR/"$m".10supers."$minc"."$maxc".split.psmcfa

# run PSMC on the genome
/share/pkg/psmc/0.6.5/bin/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o $OUTDIR/"$m".10supers."$minc"."$maxc".psmc $OUTDIR/"$m".10supers."$minc"."$maxc".psmcfa
seq 100 | xargs -i echo psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o $OUTDIR/round-{}.psmc $OUTDIR/CheMyd_draft.10supers.25.160.split.psmcfa | sh
cat $OUTDIR/${m}.10supers.25.160.psmc $OUTDIR/round-*.psmc > $OUTDIR/${m}.10supers.25.160.boots.psmc
./scripts/psmc_plot.pl -pY50000 combined $OUTDIR/${m}.10supers.25.160.boots.psmc

# plot 1y curves - to get a sense of the trend
plotdir=$OUTDIR
./scripts/psmc_plot.pl -u 1.2e-08 -g 30 "$plotdir"/"$m".10supers."$minc"."$maxc".plot "$plotdir"/"$m".10supers."$minc"."$maxc".boots.psmc
./scripts/epstopdf.pl "$plotdir"/"$m".10supers."$minc"."$maxc".boots.psmc.eps

#-----------------------------------------------------------------------------------------------------------------------
#outdir="/scratch/devel/hpawar/turtles/psmc/proc10xG/extrafilt/out"
#-rw-r--r-- 1 hpawar devel 100K Aug 25 13:50 rCheMyd1_final.10supers.18.109.psmc
#-rw-r--r-- 1 hpawar devel 100K Aug 25 14:13 rDerCor1_final.10supers.30.181.psmc

# plot 1y curves - to get a sense of the trend - both turtles on same plot - at diff y axis limits 
#plotdir=$OUTDIR

#./scripts/psmc_plot.pl -u 1.2e-08 -g 30 -Y 70  -M "CheMyd_draft.10supers.25.160.psmc,rCheMyd1_self.10supers.18.108.psmc" "$plotdir"/rCheMyd1_final_rDerCor1_final_1ycurves.plot "$outdir"/rCheMyd1_final.10supers.18.109.psmc  "$outdir"#/rDerCor1_final.10supers.30.181.psmc
#./scripts/psmc_plot.pl -u 1.2e-08 -g 30 -Y 20  -M "green,leatherback" -p "$plotdir"/rCheMyd1_final_rDerCor1_final_1ycurves_lim20.plot "$outdir"/rCheMyd1_final.10supers.18.109.psmc  "$outdir"/rDerCor1_final.10supers.30.181.psmc
#./scripts/psmc_plot.pl -u 1.2e-08 -g 30 -Y 10  -M "green,leatherback" -p "$plotdir"/rCheMyd1_final_rDerCor1_final_1ycurves_lim10.plot "$outdir"/rCheMyd1_final.10supers.18.109.psmc  "$outdir"/rDerCor1_final.10supers.30.181.psmc
