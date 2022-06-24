#########################################
### Heterozygosity when including WGR ###
#########################################

library(ggpubr)

############################
### Dermochelys coriacea ###
############################

### Whole genome ###
dc_wgr_file_list<-list.files(path = "WGR/for_genome_paper/Het/DerCor/genome_wide/")
dc_wgr_hets1<-list()
for(q in 1:length(dc_wgr_file_list)){
  df<-read.table(file=paste0("WGR/for_genome_paper/Het/DerCor/genome_wide/", dc_wgr_file_list[q]),
                 header = T)
  df$sample<-paste0("dc_",gsub("_.*","",gsub("dc_", "", dc_wgr_file_list[q])))
  df$Region<-"Whole-genome"
  dc_wgr_hets1[[q]]<-df
}
DerCor_het<-do.call(rbind, dc_wgr_hets1)
DerCor_het$het_prop<-DerCor_het$hets/DerCor_het$calls

pop_meta<-as.data.frame(levels(as.factor(DerCor_het$sample)))
pop_meta$loc<-c("US Virgin Islands","Mexico","Papua New Guinea","Costa Rica","California")
pop_meta$pop<-c("NW Atlantic","NW Atlantic","West Pacific","NW Atlantic","West Pacific")
pop_meta$ID_pop<-c("dc_NWA1","dc_NWA2","dc_WP2","dc_EP1","dc_WP1 (Reference)")
colnames(pop_meta)[1]<-"Sample_ID"

ggdensity(DerCor_het, x = "calls", y = "..density..", fill = "sample")
ggdensity(DerCor_het, x = "hets", y = "..density..", fill = "sample")
ggdensity(DerCor_het, x = "het_prop", y = "..density..", fill = "sample")
gghistogram(DerCor_het, x = "het_prop", y = "..density..", bins = 500, fill = "sample",
            xlim = c(0,0.002))

aggregate(DerCor_het$het_prop ~ DerCor_het$sample, FUN = "mean")

ggboxplot(data = DerCor_het, x = "sample", y = "het_prop", fill = "sample")

DerCor_het$loc<-pop_meta$loc[match(DerCor_het$sample,pop_meta$Sample_ID)]
DerCor_het$pop<-pop_meta$pop[match(DerCor_het$sample,pop_meta$Sample_ID)]
DerCor_het$ID_pop<-pop_meta$ID_pop[match(DerCor_het$sample,pop_meta$Sample_ID)]
DerCor_het$Spp<-"Dermochelys coriacea"

DerCor_het<-DerCor_het[order(DerCor_het$pop),]

ggdensity(DerCor_het, x = "calls", y = "..density..", fill = "ID_pop",
          xlim=c(70000,100000))
dc_boxplot_het<-ggboxplot(data = DerCor_het, x = "ID_pop", y = "het_prop", fill = "pop",
          ylim = c(0,0.00125), add = "mean", xlab = "Sample ID",
          ylab = paste0("Genome-wide heterozygosity","\n","(Proportion of heterozygotes",
                        "\n","per 100Kb window)"),
          legend.title="Population", outlier.shape = NA) +
  font("xlab", face = "bold") + font("ylab", face = "bold") + rotate_x_text(45)
ggsave(filename = "Plots/Inc_WGR/DerCor_het_w_WGR.svg", dc_boxplot_het,
       height = 16, width = 18, units = "cm")

dc_het_call_dens<-ggdensity(DerCor_het, x = "calls", y = "..density..", fill = "sample")

DerCor_genome_wide<-aggregate(DerCor_het$hets ~ DerCor_het$ID_pop, FUN = "sum")
DerCor_genome_wide[,3]<-aggregate(DerCor_het$calls ~ DerCor_het$ID_pop, FUN = "sum")[,2]
colnames(DerCor_genome_wide)<-c("Sample","Hets","Calls")
DerCor_genome_wide$het_prop<-DerCor_genome_wide$Hets/DerCor_genome_wide$Calls
mean(DerCor_genome_wide$het_prop)

### Exons ###
DerCor_exons_list<-list.files(path="WGR/for_genome_paper/Het/DerCor/exons/")
DerCor_exons_list2<-list()
for(q in 1:length(DerCor_exons_list)){
  df<-read.table(paste0("WGR/for_genome_paper/Het/DerCor/exons/",DerCor_exons_list[q]),
                 header = T)
  df$sample<-paste0("dc_",gsub("_.*","",gsub("dc_", "", DerCor_exons_list[q])))
  df$Region<-"Exons"
  DerCor_exons_list2[[q]]<-df
}
DerCor_exons<-do.call(rbind, DerCor_exons_list2)
DerCor_exons$het_prop<-DerCor_exons$hets/DerCor_exons$calls

DerCor_exons$loc<-pop_meta$loc[match(DerCor_exons$sample,pop_meta$Sample_ID)]
DerCor_exons$pop<-pop_meta$pop[match(DerCor_exons$sample,pop_meta$Sample_ID)]
DerCor_exons$ID_pop<-pop_meta$ID_pop[match(DerCor_exons$sample,pop_meta$Sample_ID)]
DerCor_exons$Spp<-"Dermochelys coriacea"

dc_exons_hets<-aggregate(DerCor_exons$hets ~ DerCor_exons$ID_pop, FUN = "sum")
dc_exons_calls<-aggregate(DerCor_exons$calls ~ DerCor_exons$ID_pop, FUN = "sum")
dc_exons<-merge(dc_exons_hets, dc_exons_calls, by = "DerCor_exons$ID_pop")
colnames(dc_exons)<-c("Sample","Hets","Calls")
dc_exons$het_prop<-dc_exons$Hets/dc_exons$Calls
mean(dc_exons$het_prop)

df<-merge(DerCor_genome_wide, dc_exons, by="Sample")
df$reduct<-df$het_prop.y/df$het_prop.x
df

### Masked genome ###
DerCor_masked_list<-list.files(path="WGR/for_genome_paper/Het/DerCor/masked/")
DerCor_masked_list2<-list()
for(q in 1:length(DerCor_masked_list)){
  df<-read.table(paste0("WGR/for_genome_paper/Het/DerCor/masked/",DerCor_masked_list[q]),
                 header = T)
  df$sample<-paste0("dc_",gsub("_.*","",gsub("dc_", "", DerCor_masked_list[q])))
  df$Region<-"masked"
  DerCor_masked_list2[[q]]<-df
}
DerCor_masked<-do.call(rbind, DerCor_masked_list2)
DerCor_masked$het_prop<-DerCor_masked$hets/DerCor_masked$calls

DerCor_masked$loc<-pop_meta$loc[match(DerCor_masked$sample,pop_meta$Sample_ID)]
DerCor_masked$pop<-pop_meta$pop[match(DerCor_masked$sample,pop_meta$Sample_ID)]
DerCor_masked$ID_pop<-pop_meta$ID_pop[match(DerCor_masked$sample,pop_meta$Sample_ID)]
DerCor_masked$Spp<-"Dermochelys coriacea"

dc_masked_hets<-aggregate(DerCor_masked$hets ~ DerCor_masked$ID_pop, FUN = "sum")
dc_masked_calls<-aggregate(DerCor_masked$calls ~ DerCor_masked$ID_pop, FUN = "sum")
dc_masked<-merge(dc_masked_hets, dc_masked_calls, by = "DerCor_masked$ID_pop")
dc_masked$het_prop<-dc_masked$`DerCor_masked$hets`/dc_masked$`DerCor_masked$calls`
mean(dc_masked$het_prop)


### Non-exons ###
DerCor_nonexons_list<-list.files(path="WGR/for_genome_paper/Het/DerCor/nonexons/")
DerCor_nonexons_list2<-list()
for(q in 1:length(DerCor_nonexons_list)){
  df<-read.table(paste0("WGR/for_genome_paper/Het/DerCor/nonexons/",DerCor_nonexons_list[q]),
                 header = T)
  df$sample<-paste0("dc_",gsub("_.*","",gsub("dc_", "", DerCor_nonexons_list[q])))
  df$Region<-"non-exons"
  DerCor_nonexons_list2[[q]]<-df
}
DerCor_nonexons<-do.call(rbind, DerCor_nonexons_list2)
DerCor_nonexons$het_prop<-DerCor_nonexons$hets/DerCor_nonexons$calls

DerCor_nonexons$loc<-pop_meta$loc[match(DerCor_nonexons$sample,pop_meta$Sample_ID)]
DerCor_nonexons$pop<-pop_meta$pop[match(DerCor_nonexons$sample,pop_meta$Sample_ID)]
DerCor_nonexons$ID_pop<-pop_meta$ID_pop[match(DerCor_nonexons$sample,pop_meta$Sample_ID)]
DerCor_nonexons$Spp<-"Dermochelys coriacea"

dc_nonexons_hets<-aggregate(DerCor_nonexons$hets ~ DerCor_nonexons$ID_pop, FUN = "sum")
dc_nonexons_calls<-aggregate(DerCor_nonexons$calls ~ DerCor_nonexons$ID_pop, FUN = "sum")
dc_nonexons<-merge(dc_nonexons_hets, dc_nonexons_calls, by = "DerCor_nonexons$ID_pop")
colnames(dc_nonexons)<-c("Sample","Hets","Calls")
dc_nonexons$het_prop<-dc_nonexons$Hets/dc_nonexons$Calls
mean(dc_nonexons$het_prop)

df<-merge(dc_nonexons, dc_exons, by="Sample")
df$reduct<-1-(df$het_prop.y/df$het_prop.x)
df


DerCor_all<-rbind(DerCor_het, DerCor_masked, DerCor_exons, DerCor_nonexons)
DerCor_all<-DerCor_all[DerCor_all$het_prop <= 0.2,]
DerCor_all<-DerCor_all[!is.na(DerCor_all$ID_pop),]

ggboxplot(DerCor_all, x = "ID_pop", y = "het_prop", outlier.shape = NA,
          ylim = c(0,0.0015), fill = "Region")

########################################################################################

######################
### Chelonia mydas ###
######################

### Whole-genome ###
cm_wgr_file_list<-list.files(path = "WGR/for_genome_paper/Het/CheMyd/whole_genome/", pattern = "")
cm_wgr_hets1<-list()
for(q in 1:length(cm_wgr_file_list)){
  df<-read.table(file=paste0("WGR/for_genome_paper/Het/CheMyd/whole_genome/", cm_wgr_file_list[q]),
                 header = T)
  df$sample<-paste0("cm_",gsub("_.*","",gsub("cm_", "", cm_wgr_file_list[q])))
  df$Region<-"Whole-genome"
  cm_wgr_hets1[[q]]<-df
}
cm_wgr_hets<-do.call(rbind, cm_wgr_hets1)
CheMyd_het<-cm_wgr_hets
CheMyd_het$het_prop<-CheMyd_het$hets/CheMyd_het$calls
CheMyd_het$Region<-"Whole-genome"


ggdensity(CheMyd_het, x = "hets", y = "..density..", fill = "sample")
ggdensity(CheMyd_het, x = "het_prop", y = "..density..", fill = "sample")
gghistogram(CheMyd_het, x = "het_prop", y = "..density..", bins = 500, fill = "sample",
            xlim = c(0,0.04))

aggregate(CheMyd_het$het_prop ~ CheMyd_het$sample, FUN = "mean")

pop_meta<-as.data.frame(levels(as.factor(CheMyd_het$sample)))
pop_meta$loc<-c("Hong Kong","Israel","South Carolina","Florida","South Carolina","Florida")
pop_meta$pop<-c("West Pacific/SE Asia","Mediterranean","NW Atlantic","NW Atlantic","NW Atlantic","NW Atlantic")
pop_meta$ID_pop<-c("cm_WP1 (Draft)","cm_MED1 (Reference)","cm_NWA1","cm_NWA2","cm_NWA3", "cm_NWA4")
colnames(pop_meta)[1]<-"Sample_ID"

CheMyd_het$loc<-pop_meta$loc[match(CheMyd_het$sample,pop_meta$Sample_ID)]
CheMyd_het$pop<-pop_meta$pop[match(CheMyd_het$sample,pop_meta$Sample_ID)]
CheMyd_het$ID_pop<-pop_meta$ID_pop[match(CheMyd_het$sample,pop_meta$Sample_ID)]
CheMyd_het$Spp<-"Chelonia mydas"

ggdensity(CheMyd_het, x = "calls", y = "..density..", fill = "ID_pop",
          xlim = c(70000,100000))
ggboxplot(data = CheMyd_het, x = "sample", y = "het_prop", fill = "sample")
ggboxplot(data = CheMyd_het, x = "ID_pop", y = "het_prop", fill = "loc",
          ylim = c(0,0.007), add = "mean", xlab = "Sample ID",
          ylab = paste0("Proportion of heterozygotes","\n","(per 100Kb window)"),
          outlier.shape = NA) +
  font("xlab", face = "bold") + font("ylab", face = "bold") + rotate_x_text(45)

CheMyd_genome_wide<-aggregate(CheMyd_het$hets ~ CheMyd_het$sample, FUN = "sum")
CheMyd_genome_wide[,3]<-aggregate(CheMyd_het$calls ~ CheMyd_het$sample, FUN = "sum")[,2]
colnames(CheMyd_genome_wide)<-c("Sample","Hets","Calls")
CheMyd_genome_wide$het_prop<-CheMyd_genome_wide$Hets/CheMyd_genome_wide$Calls
mean(CheMyd_genome_wide$het_prop)

### Exons ###
CheMyd_exons_list<-list.files(path="WGR/for_genome_paper/Het/CheMyd/exons/")
CheMyd_exons_list2<-list()
for(q in 1:length(CheMyd_exons_list)){
  df<-read.table(paste0("WGR/for_genome_paper/Het/CheMyd/exons/",CheMyd_exons_list[q]),
                 header = T)
  df$sample<-paste0("cm_",gsub("_.*","",gsub("cm_", "", CheMyd_exons_list[q])))
  df$Region<-"Exons"
  CheMyd_exons_list2[[q]]<-df
}
CheMyd_exons<-do.call(rbind, CheMyd_exons_list2)
CheMyd_exons$het_prop<-CheMyd_exons$hets/CheMyd_exons$calls


CheMyd_exons$loc<-pop_meta$loc[match(CheMyd_exons$sample,pop_meta$Sample_ID)]
CheMyd_exons$pop<-pop_meta$pop[match(CheMyd_exons$sample,pop_meta$Sample_ID)]
CheMyd_exons$ID_pop<-pop_meta$ID_pop[match(CheMyd_exons$sample,pop_meta$Sample_ID)]
CheMyd_exons$Spp<-"Chelonia mydas"


cm_exons_hets<-aggregate(CheMyd_exons$hets ~ CheMyd_exons$ID_pop, FUN = "sum")
cm_exons_calls<-aggregate(CheMyd_exons$calls ~ CheMyd_exons$ID_pop, FUN = "sum")
cm_exons<-merge(cm_exons_hets, cm_exons_calls, by = "CheMyd_exons$ID_pop")
colnames(cm_exons)<-c("Sample","Hets","Calls")
cm_exons$het_prop<-cm_exons$Hets/cm_exons$Calls
mean(cm_exons$het_prop)


### Masked genome ###
CheMyd_masked_list<-list.files(path="WGR/for_genome_paper/Het/CheMyd/masked/")
CheMyd_masked_list2<-list()
for(q in 1:length(CheMyd_masked_list)){
  df<-read.table(paste0("WGR/for_genome_paper/Het/CheMyd/masked/",CheMyd_masked_list[q]),
                 header = T)
  df$sample<-paste0("cm_",gsub("_.*","",gsub("cm_", "", CheMyd_masked_list[q])))
  df$Region<-"masked"
  CheMyd_masked_list2[[q]]<-df
}
CheMyd_masked<-do.call(rbind, CheMyd_masked_list2)
CheMyd_masked$het_prop<-CheMyd_masked$hets/CheMyd_masked$calls

CheMyd_masked$loc<-pop_meta$loc[match(CheMyd_masked$sample,pop_meta$Sample_ID)]
CheMyd_masked$pop<-pop_meta$pop[match(CheMyd_masked$sample,pop_meta$Sample_ID)]
CheMyd_masked$ID_pop<-pop_meta$ID_pop[match(CheMyd_masked$sample,pop_meta$Sample_ID)]
CheMyd_masked$Spp<-"Chelonia mydas"

cm_masked_hets<-aggregate(CheMyd_masked$hets ~ CheMyd_masked$ID_pop, FUN = "sum")
cm_masked_calls<-aggregate(CheMyd_masked$calls ~ CheMyd_masked$ID_pop, FUN = "sum")
cm_masked<-merge(cm_masked_hets, cm_masked_calls, by = "CheMyd_masked$ID_pop")
cm_masked$het_prop<-cm_masked$`CheMyd_masked$hets`/cm_masked$`CheMyd_masked$calls`
mean(cm_masked$het_prop)

### Nonexons ###
CheMyd_nonexons_list<-list.files(path="WGR/for_genome_paper/Het/CheMyd/nonexons/")
CheMyd_nonexons_list2<-list()
for(q in 1:length(CheMyd_nonexons_list)){
  df<-read.table(paste0("WGR/for_genome_paper/Het/CheMyd/nonexons/",CheMyd_nonexons_list[q]),
                 header = T)
  df$sample<-paste0("cm_",gsub("_.*","",gsub("cm_", "", CheMyd_nonexons_list[q])))
  df$Region<-"non-exons"
  CheMyd_nonexons_list2[[q]]<-df
}
CheMyd_nonexons<-do.call(rbind, CheMyd_nonexons_list2)
CheMyd_nonexons$het_prop<-CheMyd_nonexons$hets/CheMyd_nonexons$calls

CheMyd_nonexons$loc<-pop_meta$loc[match(CheMyd_nonexons$sample,pop_meta$Sample_ID)]
CheMyd_nonexons$pop<-pop_meta$pop[match(CheMyd_nonexons$sample,pop_meta$Sample_ID)]
CheMyd_nonexons$ID_pop<-pop_meta$ID_pop[match(CheMyd_nonexons$sample,pop_meta$Sample_ID)]
CheMyd_nonexons$Spp<-"Chelonia mydas"

cm_nonexons_hets<-aggregate(CheMyd_nonexons$hets ~ CheMyd_nonexons$ID_pop, FUN = "sum")
cm_nonexons_calls<-aggregate(CheMyd_nonexons$calls ~ CheMyd_nonexons$ID_pop, FUN = "sum")
cm_nonexons<-merge(cm_nonexons_hets, cm_nonexons_calls, by = "CheMyd_nonexons$ID_pop")
colnames(cm_nonexons)<-c("Sample","Hets","Calls")
cm_nonexons$het_prop<-cm_nonexons$Hets/cm_nonexons$Calls
mean(cm_nonexons$het_prop)

df<-merge(cm_nonexons, cm_exons, by="Sample")
df$reduct<-1-(df$het_prop.y/df$het_prop.x)
df
mean(df$reduct)

CheMyd_all<-rbind(CheMyd_het, CheMyd_masked, CheMyd_exons, CheMyd_nonexons)
CheMyd_all<-CheMyd_all[CheMyd_all$het_prop <= 0.2,]
CheMyd_all<-CheMyd_all[!is.na(CheMyd_all$ID_pop),]

ggboxplot(CheMyd_all, x = "ID_pop", y = "het_prop", outlier.shape = NA,
          ylim = c(0,0.0075), fill = "Region")

### Combined ###
all_het<-rbind(DerCor_all, CheMyd_all)
all_het<-all_het[order(all_het$Region),]
A<-ggboxplot(all_het, x = "Region", y = "het_prop", fill = "ID_pop",
          outlier.shape = NA, ylim = c(0,0.006), title = "A")
B<-ggboxplot(all_het, x = "ID_pop", y = "het_prop", fill = "Region",
             outlier.shape = NA, ylim = c(0,0.006), title = "B") +  rotate_x_text(45)
C<-ggarrange(A, B, ncol = 1, nrow = 2)
ggsave(filename = "Plots/Inc_WGR/Sample_het_plot.png", C, height = 21, width = 21, units = "cm")
ggsave(filename = "Plots/Inc_WGR/heterozygosity.svg", B, height = 18, width = 21, units = "cm")

exons<-as.data.frame(rbind(DerCor_exons, CheMyd_exons))
exons<-exons[exons$calls >= 100,]
meanex<-mean(exons$het_prop)
exons_high<-exons[]

#######################

zerowins<-all_het[all_het$Region == "masked" &
                    all_het$het_prop == 0,]
zerowin_t<-as.data.frame(table(zerowins$ID_pop))
zerowin_t$Spp<-c(rep("Chelonia", 6), rep("Dermochelys", 5))
t.test(zerowin_t$Freq ~ zerowin_t$Spp)
zeroplot<-ggboxplot(zerowin_t, x = "Spp", y = "Freq", fill = "Spp")
ggsave(filename = "Plots/Inc_WGR/zero_het_box.svg", zeroplot,
       height = 18, width = 18, units = "cm")


DerCor_genome_wide$Region<-"Genome-wide"
dc_masked$Region<-"Masked"
colnames(dc_masked)<-c("Sample","Hets","Calls","het_prop","Region")
dc_exons$Region<-"Exons"
dc_nonexons$Region<-"Nonexons"
dc<-rbind(DerCor_genome_wide, dc_masked, dc_exons, dc_nonexons)
dc$Spp<-"Dermochelys"

CheMyd_genome_wide$Region<-"Genome-wide"
cm_masked$Region<-"Masked"
colnames(cm_masked)<-c("Sample","Hets","Calls","het_prop","Region")
cm_exons$Region<-"Exons"
cm_nonexons$Region<-"Nonexons"
cm<-rbind(CheMyd_genome_wide, cm_masked, cm_exons, cm_nonexons)
cm$Spp<-"Chelonia"

all_regions<-rbind(dc,cm)

masked<-all_regions[all_regions$Region == "Masked",]
t.test(masked$het_prop ~ masked$Spp)

exons<-all_regions[all_regions$Region == "Exons",]
t.test(exons$het_prop ~ exons$Spp)

dcexno<-dc[dc$Region == "Exons" |
             dc$Region == "Nonexons",]
t.test(dcexno$het_prop ~ dcexno$Region, paired = T)

cmexno<-cm[cm$Region == "Exons" |
             cm$Region == "Nonexons",]
t.test(cmexno$het_prop ~ cmexno$Region, paired = T)
