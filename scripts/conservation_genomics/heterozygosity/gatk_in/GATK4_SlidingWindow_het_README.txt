GATK4_SlidingWindow_heterozygosity

Scripts to detect variants in a genome assembly and count them across non-overlapping windows of specified size (e.g., 1MB).

These are adapted from methods originally described in Robinson, J.A., Raikkonen, J., Vucetich, L.M., Vucetich, J.A., Peterson, R.O., Lohmueller, K.E., Wayne, R.K., 2019. Genomic signatures of extensive inbreeding in Isle Royale wolves, a population on the threshold of extinction. Science Advances 5, eaau0757.

SCRIPTS:
gatk4_pipeline_winhet_300820.sh
slidingWindowHet_280820.py
filterVCF_300820.py

Required files:
Alignment of sequence reads to a reference sequence (.bam file). This script assumes one sample per bam file, so each sample alignment to a reference should be run through the pipeline separately (this differs from J. Robinson's original pipeline, which handled multiple merged alignments from different samples). 

Reference sequence (.fasta file, indexed with samtools index)

scaffolds list (text file listing scaffolds in the reference fasta file). See below for how to generate scaffold list

chrom_lengths.txt (text file listing scaffold names and their lengths (tab separated)). See below for how to generate chromosome lengths list


Required values:
WINSIZE=1000000
STEPSIZE=1000000
MINDEPTH (integer; minimum depth for calling genotypes)
MAXDEPTH (integer; maximum depth for calling genotypes)
	# The minimum depth is typically 1/3 average depth of coverage; max is 2x average depth of coverage. 
##########################
#### Scaffolds list file ####
Generate the scaffold list for the reference genome (from the fasta file)

grep "^>" GCF_008692025.1_mPhoSin1.pri_genomic.fna > GCF_008692025.1_mPhoSin1.pri_genomic.fna.scaffolds.list

# This list has the ">" at the beginning of each line, so need to remove those in a text editor, then save and copy back to folder with the reference fasta file. 
# There's probably an easier way to do this with awk (see below), but I haven't figured it out.

#########################
Generate chromosome lengths from the fasta index file (from samtools index):
awk '{print $1}' filename.fasta.fai > filename.chrom_lengths.txt



####################################################################
Running script
Set the paths to files and scripts in the shell script gatk4_pipeline_winhet_280820.sh
Determine which scaffolds to scan for SNPs and set them as the array numbers (e.g., in the sbatch header: #SBATCH --array=1-20,22,24

##########################
OUTPUT files
1) The filtered VCF file (*${NUM}.Filter.vcf.gz) contains the final genotype information for every site in the alignment. VCF files are generated for each part of the GATK pipeline, but since the last step (filter.vcf) is the cumulative effect of the pipeline steps, it's the only one stored in the output directory "OUTDIR". The others are stored in TEMPDIR, which also speeds the script if writing those files is on the scratch directory. 

2) Sliding window genotype call and heterozygote counts for each scaffold are stored in the scaffold-specific files *${NUM}.Filter.vcf.gz_het_1000000win_1000000step.txt

##########################
Summarizing and plotting data
The R script "winHet_plot_het_data_010920.R". 
The R script loads the files based on filenames based on the window size (e.g., 10000000) and "step.txt", and filters based on the specified scaffold number (based on order in scaffolds list). 
OUTPUT files = 4 pdf files



