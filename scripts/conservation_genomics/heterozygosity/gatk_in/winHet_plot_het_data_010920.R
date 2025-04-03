library(plyr)
rm(list = ls())

# setwd("~/Dropbox/Research/Vaquita/analysis/single_genomes/winHet")
#setwd("~/Documents/Mol_Ecol_Lab/_Projects/ blue_whales/Bmus_Ref_Genome_Morgridge_VGP/Bmus_Illumina_Yury/Bmus_geno_ROH_scripts/ROH_plots")
getwd()
######################################################################
ID<-"Psin" # Project ID (added to output file names)
Scaffolds<-c(1:7,9:22) # includes scaffolds ≥10 Mb in length;
# exclude X (file 34, "scaffold_4_arrow_ctg1".
date<-format(Sys.Date(), "%d.%m.%y")

# z017954, 1 Mb windows, exclude X chromosome (scaffold_4_arrow) and scaffolds <1 Mb in length
winhetfiles=list.files(pattern="_1000000step.txt")
winhetfiles=winhetfiles[Scaffolds]     # [c(1:3,5:29)]
# remove NA's 
winhetfiles<-winhetfiles[!is.na(winhetfiles)]
allhet=ldply(winhetfiles, read.table, header=TRUE, sep="\t")

# overall heterozygosity (hets/calls across all scaffolds)
sumHets<-sum(allhet$hets)
sumCalls<-sum(allhet$calls)
Het_1Mb<-round(sumHets/sumCalls,7)
allhet$pct_het<-allhet$hets/allhet$calls
sd_het <- sd(allhet$pct_het)

# calculate minimum and maximum heterozygosity from all 1Mb intervals
minHet_1Mb<-min(allhet$hets, na.rm=T)
maxhet_count<-max(allhet$hets)
maxcalls_ID<-which.max(allhet$hets)
Calls_maxhet<-allhet[maxcalls_ID,4]
Max_het_1MB<-maxhet_count/Calls_maxhet

# total number of heterozygous sites
tot_hets <- sum(allhet$hets)

# Count number of windows with heterozygosity = 0
count_nohet<-sum(allhet$hets == 0)
# percent of windows with no hets
percent_nohet<-count_nohet/nrow(allhet)
# count number of windows with heterozygosity < 0.0001 (0.1/kb = first bar in windows plot)
count_het_0.1_per_kb<-sum(allhet$pct_het < 0.0001)

# print out het count summary
a=paste0("Min. heterozygote count/Mb = ",minHet_1Mb)
b=paste0("Max. heterozygote count/Mb = ",maxhet_count)
c=paste0("Total number of heterozygote site = ", tot_hets)
d=paste0("Number of 1Mb windows with het=0 ",count_nohet)
e=paste0("Proportion of 1Mb windows with het=0 ",percent_nohet)
f=paste0("Max heterozygosity = ", Max_het_1MB)
g=paste0("Average heterozygosity = ", Het_1Mb)
h=paste0("Standard deviation of het/window = ",sd_het)

pdf(paste0(ID, "_het_summ_stats.pdf"))
plot(NA, xlim=c(0,8), ylim=c(0,8), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,7,a, pos=4)
text(1,6,b, pos=4)
text(1,5,f, pos=4)
text(1,4,c, pos=4)
text(1,3,d, pos=4)
text(1,1,e, pos=4)
text(1,8,g, pos=4)
text(1,2,h, pos=4)
points(rep(1,8),1:8, pch=15)
dev.off()


# Get boundaries of chromosomes for plotting
pos=as.numeric(rownames(unique(data.frame(allhet$chrom)[1])))
pos=append(pos,length(allhet$chrom))
numpos=NULL
for (i in 1:length(pos)-1){numpos[i]=(pos[i]+pos[i+1])/2}

# Set plot colors (2 colors alternating between chromosomes)
mycols=NULL
for (i in (seq(1,length(numpos), by=2))){mycols[i]="#2171b5"}
for (i in (seq(2,length(numpos), by=2))){mycols[i]="#6baed6"}

# Barplot of heterozygosity in windows across the genome
pdf(paste0(ID,"_winhet_barplot_1Mbwin_",date,".pdf"), width=8, height=4, pointsize=10)
par(mar=c(8,5,2,1))
b=barplot(1000*allhet$hets/allhet$calls, ylim=c(0,6),
          border=mycols[as.factor(allhet$chrom)], col=mycols[as.factor(allhet$chrom)], 
          ylab="Het / kb", , cex.lab=1.25, main=paste(ID,"(Het=",Het_1Mb*100,"%)"))
axis(side=1, at=b[pos], labels=F)
axis(side=1, at=b[numpos], tick=F, labels=as.character(unique(allhet$chrom)), las=3, 
     line=-.25, cex.axis=.8)
dev.off()

# Barplot of heterozygosity in windows across the genome, without main label, scaffold labels
pdf(paste0(ID,"_winhet_barplot_1Mbwin_NoTitle",date,".pdf"), width=8, height=4, pointsize=14)
par(mar=c(8,5,2,1))
b=barplot(1000*allhet$hets/allhet$calls, ylim=c(0,6),
          border=mycols[as.factor(allhet$chrom)], col=mycols[as.factor(allhet$chrom)], 
          ylab="Het / kb", cex.lab=1.25)
axis(side=1, at=b[pos], labels=F)
axis(side=1, at=b[numpos], tick=F, labels=as.character(unique(allhet$chrom)), las=3, 
     line=-.25, cex.axis=0.8)
dev.off()

# Histogram of per-window heterozygosity
pdf(paste0(ID,"_winhet_hist_1Mbwin_",date,".pdf"), width=4, height=4, pointsize=14)
par(mar=c(5,5,2,1))
h=hist(1000*allhet$hets/allhet$calls, breaks=seq(0,6, by=0.1), ylim=c(0,200), 
       border="#2171b5", col="#2171b5", ylab="# of windows", xlab="Het / kb", cex.lab=1.5,
       main=paste(ID,"(Het=",Het_1Mb*100,"%)"))
dev.off()

