Scripts in this directory are directly relevant to the conservation genomics component of the genome analysis. This includes:
  1. Alignment of the processed 10X reads (see pre-processing steps) to the reference genomes (BWA-mem).
  2. Removal of PCR duplicates (Picard Tools).
  3. Addition of read group headers (Picard Tools).
  4. Re-alignment around indels (GATK)*.
  5. Heterozygosity estimates, specifically the use of the GATK pipeline to estimate heterozygosity for:
    a. genome-wide (GATK, ANGSD).
    b. repeat and low-complexity masked genome (see transposable element steps for generation of the mask file) (GATK).
    c. only regions identified as exons in the annotation of the genomes (GATK).
    d. non-exon regions, essentially the inverse of 'c' above - only regions where no exons are present (GATK).
  6. Analysis of runs of homozygosity (ROHs)
    a. Trimming, alignment, duplicate removal, and indel realignment of low/medium-coverage whole-genome resequence data (Trimmomatic, BWA-mem, Picard Tools, GATK).
    b. Generation of a SNP-list for use with PLINK (ANGSD).
    c. ROH analysis using the SNP-based PLINK method (PLINK).
   7. Analysis of accumulation of deleterious alleles and genetic load.
   8. Exploration of genes of interest for conservation.
    a. Genes related to temperature-dependent sex-determination.
    b. Genes associated the immune system.
    
    
* Note, for Step 4: require GATK version prior to 4.0, where indel realignment was incoporated, rather than a stand alone function.
