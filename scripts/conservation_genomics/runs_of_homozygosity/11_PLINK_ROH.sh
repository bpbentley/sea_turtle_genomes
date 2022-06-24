#!/bin/bash
#BSUB -J CheMyd_ROH_PLINK
#BSUB -o ./logs/CheMyd_ROH_PLINK.log
#BSUB -e ./logs/CheMyd_ROH_PLINK.err
#BSUB -n 20
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=6000]
#BSUB -W 04:00
#BSUB -q short

####################
### Load modules ###
####################

module load plink/1.90b6.9

######################
### Set file paths ###
######################
INDIR=/project/uma_lisa_komoroske/Blair/WGR_for_genome/ANGSD/CheMyd
PLINKDIR=/project/uma_lisa_komoroske/Blair/WGR_for_genome/ROH_PLINK/CheMyd


#####################################################
### Run PLINK ROH estimation for each chromosome: ###
#####################################################

for p in {1..28}; do
plink --tfile $INDIR/CheMyd_SUPER_${p} --homozyg-snp 20 --homozyg-kb 1 \
 --homozyg-window-snp 20 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.01 \
 --allow-extra-chr --out $PLINKDIR/CheMyd_SUPER_${p}
done

##########################
### Testing Paramaters ###
##########################

#for q in {0..100..10}; do
#for p in {1..28}; do
#plink --tfile $INDIR/'subset_DerCor_SUPER_'$p --homozyg-snp 20 --homozyg-kb 50 \
# --homozyg-window-snp 20 --homozyg-window-het 3 --homozyg-window-missing ${q} --homozyg-window-threshold 0.01 \
# --allow-extra-chr --out $PLINKDIR/test_MISS/MISS_${q}/subset_DerCor_${p}_MISS_${q}
#done
#done

# --homozyg-density 20 --homozyg-gap 1000

#for p in {1..28}; do
#for q in {10..200..10}; do
#plink --tfile $INDIR/'WGR_DerCorTop5_SUPER_'${p} --homozyg-snp ${q} --homozyg-kb 1 \
# --homozyg-window-snp ${q} --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.01 \
# --allow-extra-chr --out $PLINKDIR/test_SNP/SNP_${q}/DerCorTop5_SUPER_${p}_SNP_${q}
#done
#done