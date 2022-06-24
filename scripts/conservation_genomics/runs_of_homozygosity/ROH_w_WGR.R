#################################
### ROH Analysis after review ###
#################################

library(ggpubr)
library(mclust)

### Determine size categories for ROHs based on Pemberton et al. 2012
test1<-Mclust(all$KB, G = 3)
summary(test1, parameters = T)
plot.Mclust(test1, what = "classification")
ROH_classes<-as.data.frame(cbind(test1$data, test1$classification))
aggregate(ROH_classes$V1 ~ ROH_classes$V2, FUN = "max")
aggregate(ROH_classes$V1 ~ ROH_classes$V2, FUN = "min")

# Individual species?
dcHC<-Mclust(dc_df$KB, G = 3)
summary(test1, parameters = T)
dcROH_classes<-as.data.frame(cbind(dcHC$data, dcHC$classification))
aggregate(dcROH_classes$V1 ~ dcROH_classes$V2, FUN = "max")
aggregate(dcROH_classes$V1 ~ dcROH_classes$V2, FUN = "min")

cmHC<-Mclust(cm_df$KB, G = 3)
summary(test1, parameters = T)
cmROH_classes<-as.data.frame(cbind(cmHC$data, cmHC$classification))
aggregate(cmROH_classes$V1 ~ cmROH_classes$V2, FUN = "max")
aggregate(cmROH_classes$V1 ~ cmROH_classes$V2, FUN = "min")

## DerCor
dc_samples<-c("dc_101533","dc_11171","dc_20292","dc_33126","rDerCor1")
num_samples<-c(seq(1:5))
dc_samples<-as.data.frame(cbind(num_samples,dc_samples))
dc_samples$pop_ID<-c("dc_NWA1","dc_EP1","dc_WP2","dc_NWA2","dc_WP1 (Reference)")
dc_list<-list.files(path="WGR/for_genome_paper/ROH/DerCor")
dc_array<-list()
for(q in 1:28){
  df<-read.table(file = paste0("WGR/for_genome_paper/ROH/DerCor/",dc_list[q]),
                 header = T)
  dc_array[[q]]<-df
}
dc_df<-do.call(rbind, dc_array)
dc_df$SID<-dc_samples$pop_ID[match(dc_df$FID, dc_samples$num_samples)]
dc_df<-dc_df[dc_df$KB >= 50,]

pop_meta<-as.data.frame(levels(as.factor(dc_df$SID)))
pop_meta$loc<-c("US Virgin Islands","Mexico","Papua New Guinea","Costa Rica","California")
pop_meta$pop<-c("NW Atlantic","NW Atlantic","West Pacific","NW Atlantic","West Pacific")
pop_meta$ID_pop<-c("dc_NWA1","dc_EP1","dc_WP2","dc_NWA2","dc_WP1 (Reference)")

ggboxplot(dc_df, x = "SID", y = "KB")
dc_df$logKB<-log10(dc_df$KB)
ggboxplot(dc_df, x = "SID", y = "logKB", fill = "SID",
          xlab = "Sample ID", ylab = "log(Length of ROH)", legend = "none") +
  font("xlab", face = "bold") + font("ylab", face = "bold")
ggviolin(dc_df, x = "SID", y = "logKB", fill = "SID")

aggregate(dc_df$KB ~ dc_df$FID, FUN = "mean")
aggregate(dc_df$KB ~ dc_df$FID, FUN = "max")

dc_small<-dc_df[dc_df$KB <= 500,]
dc_small$Category<-"Short"
dc_small_cat<-aggregate(dc_small$KB ~ dc_small$SID, FUN = "sum")
colnames(dc_small_cat)<-c("Sample_ID", "sumKB")
dc_small_cat$Category<-"Short"

dc_medium<-dc_df[dc_df$KB > 500 & dc_df$KB <= 1000,]
dc_medium$Category<-"Medium"
dc_medium_cat<-aggregate(dc_medium$KB ~ dc_medium$SID, FUN = "sum")
colnames(dc_medium_cat)<-c("Sample_ID", "sumKB")
dc_medium_cat$Category<-"Medium"

dc_long<-dc_df[dc_df$KB > 1000,]
dc_long$Category<-"Long"
dc_long_cat<-aggregate(dc_long$KB ~ dc_long$SID, FUN = "sum")
colnames(dc_long_cat)<-c("Sample_ID", "sumKB")
dc_long_cat$Category<-"Long"

dc_df<-as.data.frame(rbind(dc_small, dc_medium, dc_long))

dc_cat<-rbind(dc_small_cat,dc_medium_cat,dc_long_cat)
dc_cat$sumMB<-dc_cat$sumKB/1000

ggbarplot(dc_cat, x = "Sample_ID", y = "sumMB", fill = "Category",
          xlab = "Sample ID", ylab = "Aggregate length of ROHs (MB)",
          legend.title="ROH Category", orientation = "horiz") +
  font("xlab", face = "bold") + font("ylab", face = "bold")

dc_sum<-aggregate(dc_df$KB ~ dc_df$SID, FUN = "sum")
colnames(dc_sum)<-c("Sample_ID","Agg_length")
dc_sum$Count<-table(dc_df$SID)


ggdensity(data = dc_df, x = "logKB", y = "..density..", fill = "SID",
          facet.by = "SID")
table(dc_df$SID)
aggregate(dc_df$KB ~ dc_df$SID, FUN = "sum")

summary(aov(dc_df$KB ~ dc_df$Category))

## CheMyd
cm_samples<-c("Flower","Yucca","Poppy","27-2017-Cm","CheMyd_draft","rCheMyd1")
num_samples<-c(seq(1:6))
cm_samples<-as.data.frame(cbind(num_samples,cm_samples))
cm_samples$pop_ID<-c("cm_NWA1","cm_NWA2", "cm_NWA3", "cm_NWA4", "cm_WP1 (Draft)", "cm_MED1 (Reference)")
cm_list<-list.files(path="WGR/for_genome_paper/ROH/CheMyd")
cm_array<-list()
for(q in 1:28){
  df<-read.table(file = paste0("WGR/for_genome_paper/ROH/CheMyd/",cm_list[q]),
                 header = T)
  cm_array[[q]]<-df
}
cm_df<-do.call(rbind, cm_array)
cm_df$SID<-cm_samples$pop_ID[match(cm_df$FID, cm_samples$num_samples)]
cm_df<-cm_df[cm_df$KB >= 50,]

ggboxplot(cm_df, x = "SID", y = "KB")
cm_df$logKB<-log10(cm_df$KB)
ggboxplot(cm_df, x = "SID", y = "logKB", fill = "SID",
          xlab = "Sample ID", ylab = "log(Length of ROH)", legend = "none") +
  font("xlab", face = "bold") + font("ylab", face = "bold") + rotate_x_text(45)
ggviolin(cm_df, x = "SID", y = "logKB", fill = "SID")

aggregate(cm_df$KB ~ cm_df$SID, FUN = "mean")
aggregate(cm_df$KB ~ cm_df$SID, FUN = "max")
aggregate(cm_df$KB ~ cm_df$SID, FUN = "min")

cm_small<-cm_df[cm_df$KB <= 500,]
cm_small$Category<-"Short"
cm_small_cat<-aggregate(cm_small$KB ~ cm_small$SID, FUN = "sum")
colnames(cm_small_cat)<-c("Sample_ID", "sumKB")
cm_small_cat$Category<-"Short"

cm_medium<-cm_df[cm_df$KB > 500 & cm_df$KB <= 1000,]
cm_medium$Category<-"Medium"
cm_medium_cat<-aggregate(cm_medium$KB ~ cm_medium$SID, FUN = "sum")
colnames(cm_medium_cat)<-c("Sample_ID", "sumKB")
cm_medium_cat$Category<-"Medium"

cm_long<-cm_df[cm_df$KB > 1000,]
cm_long$Category<-"Long"
cm_long_cat<-aggregate(cm_long$KB ~ cm_long$SID, FUN = "sum")
colnames(cm_long_cat)<-c("Sample_ID", "sumKB")
cm_long_cat$Category<-"Long"

cm_df<-as.data.frame(rbind(cm_small,cm_medium,cm_long))

cm_cat<-rbind(cm_small_cat,cm_medium_cat,cm_long_cat)
cm_cat$sumMB<-cm_cat$sumKB/1000
ggbarplot(cm_cat, x = "Sample_ID", y = "sumMB", fill = "Category",
          xlab = "Sample ID", ylab = "Aggregate length of ROHs (MB)",
          legend.title="ROH Category", orientation = "horiz") +
  font("xlab", face = "bold") + font("ylab", face = "bold")

cm_sum<-aggregate(cm_df$KB ~ cm_df$SID, FUN = "sum")
colnames(cm_sum)<-c("Sample_ID","Agg_length")
cm_sum$Count<-table(cm_df$SID)

ggdensity(data = cm_df, x = "logKB", y = "..density..", fill = "SID",
          facet.by = "SID")
table(cm_df$SID)
aggregate(cm_df$KB ~ cm_df$SID, FUN = "sum")

## Both
dc_df$Spp<-"DerCor"
cm_df$Spp<-"CheMyd"
all<-rbind(cm_df, dc_df)
dc_df$Chrom<-as.numeric(gsub("SUPER_", "", dc_df$CHR))
dc_df<-dc_df[order(dc_df$Chrom),]
dc_pal<-get_palette(palette = "Blues", 5)

cm_df$Chrom<-as.numeric(gsub("SUPER_", "", cm_df$CHR))
cm_df<-cm_df[order(cm_df$Chrom),]
cm_pal<-get_palette(palette = "Reds", 6)

dc_ROHs<-ggboxplot(dc_df, x = "CHR", y = "logKB", fill = "SID",
          xlab = "Chromosome", ylab = "Length of ROH (logKB)",
          legend.title = "ID", palette = dc_pal, title = "(a)") + rotate_x_text(45)

cm_ROHS<-ggboxplot(cm_df, x = "CHR", y = "logKB", fill = "SID",
          xlab = "Chromosome", ylab = "Length of ROH (logKB)",
          legend.title = "ID", palette = cm_pal, title = "(b)") + rotate_x_text(45)
ROHs_plot<-ggarrange(dc_ROHs, cm_ROHS, ncol = 1, nrow = 2)  
ggsave(filename = "Plots/Inc_WGR/ROHS_across_genome.png", height = 21, width = 29, units = "cm")


gghistogram(all, x = "logKB", y = "..density..", bins = 100, fill = "Spp")
all$lnKB<-log(all$KB)
gghistogram(all, x = "lnKB", y = "..density..", bins = 100, fill = "Spp")

t.test(all$KB ~ all$Spp)
all_small<-all[all$Category == "Short",]
summary(aov(all_small$KB ~ all_small$Spp))

all_cat<-as.data.frame(rbind(dc_cat, cm_cat))

all_cat<-all_cat[order(all_cat$Sample_ID, decreasing = T),]
fig4c<-ggbarplot(all_cat, x = "Sample_ID", y = "sumMB", fill = "Category",
          xlab = "Sample ID", ylab = "Aggregate length of ROHs (MB)",
          legend.title="ROH Category", orientation = "horiz") +
  font("xlab", face = "bold") + font("ylab", face = "bold")
ggsave(file="Plots/Inc_WGR/fig4c.svg", fig4c, height = 18, width = 18, units = "cm")

number<-as.data.frame(table(all$SID))
colnames(number)<-c("V1","V2")
number$V3<-"Number of ROHs"
length<-aggregate(all$KB ~ all$SID, FUN = "sum")
colnames(length)<-c("V1","V2")
length$V3<-"Total length of ROHs"
both<-merge(number, length, by = "V1")
both$spp<-c(rep("Chelonia",6), rep("Dermochelys",5))
ggscatter(both, x = "V2.y", y = "V2.x", col = "spp",
          ylab = "Total number of ROHS (>50Kb)",
          xlab = "Total aggregate length of ROHs (kb)",
          legend.title = "Species") +
  font("xlab", face = "bold") + font("ylab", face = "bold")

### CheMyd aligned to Draft
cmd_samples<-c("CheMyd_draft","rCheMyd1","Flower","Yucca","Poppy","27-2017-Cm")
num_samples<-c(seq(1:6))
cmd_samples<-as.data.frame(cbind(num_samples,cmd_samples))
cmd_samples$pop_ID<-c("cm_WP1 (Draft)","cm_MED1 (Reference)","cm_NWA1","cm_NWA2", "cm_NWA3", "cm_NWA4")
cmd_list<-list.files(path="WGR/for_genome_paper/ROH/draft_CheMyd")
cmd_array<-list()
for(q in 1:518){
  df<-read.table(file = paste0("WGR/for_genome_paper/ROH/draft_CheMyd/",cmd_list[q]),
                 header = T)
  cmd_array[[q]]<-df
}
cmd_df<-do.call(rbind, cmd_array)
cmd_df$SID<-cmd_samples$pop_ID[match(cmd_df$FID, cmd_samples$num_samples)]
cmd_df<-cmd_df[cmd_df$KB >= 50,]

ggboxplot(cmd_df, x = "SID", y = "KB")
cmd_df$logKB<-log10(cmd_df$KB)
ggboxplot(cmd_df, x = "SID", y = "logKB", fill = "SID",
          xlab = "Sample ID", ylab = "log(Length of ROH)", legend = "none") +
  font("xlab", face = "bold") + font("ylab", face = "bold") + rotate_x_text(45)
ggviolin(cmd_df, x = "SID", y = "logKB", fill = "SID")

aggregate(cmd_df$KB ~ cmd_df$SID, FUN = "mean")
aggregate(cmd_df$KB ~ cmd_df$SID, FUN = "max")
aggregate(cmd_df$KB ~ cmd_df$SID, FUN = "min")

cmd_small<-cmd_df[cmd_df$KB <= 500,]
cmd_small_cat<-aggregate(cmd_small$KB ~ cmd_small$SID, FUN = "sum")
colnames(cmd_small_cat)<-c("Sample_ID", "sumKB")
cmd_small_cat$Category<-"Short"

cmd_medium<-cmd_df[cmd_df$KB > 500 & cmd_df$KB <= 1000,]
cmd_medium_cat<-aggregate(cmd_medium$KB ~ cmd_medium$SID, FUN = "sum")
colnames(cmd_medium_cat)<-c("Sample_ID", "sumKB")
cmd_medium_cat$Category<-"Medium"

cmd_long<-cmd_df[cmd_df$KB > 1000,]
cmd_long_cat<-aggregate(cmd_long$KB ~ cmd_long$SID, FUN = "sum")
colnames(cmd_long_cat)<-c("Sample_ID", "sumKB")
cmd_long_cat$Category<-"Long"


cmd_cat<-rbind(cmd_small_cat,cmd_medium_cat)
cmd_cat$sumMB<-cmd_cat$sumKB/1000
cmd_cat$Reference<-"Draft"
ggbarplot(cmd_cat, x = "Sample_ID", y = "sumMB", fill = "Category",
          xlab = "Sample ID", ylab = "Aggregate length of ROHs (MB)",
          legend.title="ROH Category", orientation = "horiz") +
  font("xlab", face = "bold") + font("ylab", face = "bold")

cmd_sum<-aggregate(cmd_df$KB ~ cmd_df$SID, FUN = "sum")
colnames(cmd_sum)<-c("Sample_ID","Agg_length")
cmd_sum$Count<-table(cmd_df$SID)

ggdensity(data = cmd_df, x = "logKB", y = "..density..", fill = "SID",
          facet.by = "SID")
table(cmd_df$SID)
aggregate(cmd_df$KB ~ cmd_df$SID, FUN = "sum")

## CheMyd aligned to DNAZoo version

cmzoo_samples<-c("CheMyd_draft","rCheMyd1","Flower","Yucca","Poppy","27-2017-Cm")
num_samples<-c(seq(1:6))
cmzoo_samples<-as.data.frame(cbind(num_samples,cmzoo_samples))
cmzoo_samples$pop_ID<-c("cm_WP1 (Draft)","cm_MED1 (Reference)","cm_NWA1","cm_NWA2", "cm_NWA3", "cm_NWA4")
cmzoo_list<-list.files(path="WGR/for_genome_paper/ROH/DNAzoo_CheMyd/")
cmzoo_array<-list()
for(q in 1:28){
  df<-read.table(file = paste0("WGR/for_genome_paper/ROH/DNAzoo_CheMyd/",cmzoo_list[q]),
                 header = T)
  cmzoo_array[[q]]<-df
}
cmzoo_df<-do.call(rbind, cmzoo_array)
cmzoo_df$SID<-cmzoo_samples$pop_ID[match(cmzoo_df$FID, cmzoo_samples$num_samples)]
cmzoo_df<-cmzoo_df[cmzoo_df$KB >= 50,]

ggboxplot(cmzoo_df, x = "SID", y = "KB")
cmzoo_df$logKB<-log10(cmzoo_df$KB)
ggboxplot(cmzoo_df, x = "SID", y = "logKB", fill = "SID",
          xlab = "Sample ID", ylab = "log(Length of ROH)", legend = "none") +
  font("xlab", face = "bold") + font("ylab", face = "bold") + rotate_x_text(45)
ggviolin(cmzoo_df, x = "SID", y = "logKB", fill = "SID")

aggregate(cmzoo_df$KB ~ cmzoo_df$SID, FUN = "mean")
aggregate(cmzoo_df$KB ~ cmzoo_df$SID, FUN = "max")
aggregate(cmzoo_df$KB ~ cmzoo_df$SID, FUN = "min")

cmzoo_small<-cmzoo_df[cmzoo_df$KB <= 500,]
cmzoo_small_cat<-aggregate(cmzoo_small$KB ~ cmzoo_small$SID, FUN = "sum")
colnames(cmzoo_small_cat)<-c("Sample_ID", "sumKB")
cmzoo_small_cat$Category<-"Short"

cmzoo_medium<-cmzoo_df[cmzoo_df$KB > 500 & cmzoo_df$KB <= 1000,]
cmzoo_medium_cat<-aggregate(cmzoo_medium$KB ~ cmzoo_medium$SID, FUN = "sum")
colnames(cmzoo_medium_cat)<-c("Sample_ID", "sumKB")
cmzoo_medium_cat$Category<-"Medium"

cmzoo_long<-cmzoo_df[cmzoo_df$KB > 1000,]
cmzoo_long_cat<-aggregate(cmzoo_long$KB ~ cmzoo_long$SID, FUN = "sum")
colnames(cmzoo_long_cat)<-c("Sample_ID", "sumKB")
cmzoo_long_cat$Category<-"Long"


cmzoo_cat<-rbind(cmzoo_small_cat,cmzoo_medium_cat)
cmzoo_cat$sumMB<-cmzoo_cat$sumKB/1000
cmzoo_cat$Reference<-"DNAZoo"
ggbarplot(cmzoo_cat, x = "Sample_ID", y = "sumMB", fill = "Category",
          xlab = "Sample ID", ylab = "Aggregate length of ROHs (MB)",
          legend.title="ROH Category", orientation = "horiz") +
  font("xlab", face = "bold") + font("ylab", face = "bold")

cmzoo_sum<-aggregate(cmzoo_df$KB ~ cmzoo_df$SID, FUN = "sum")
colnames(cmzoo_sum)<-c("Sample_ID","Agg_length")
cmzoo_sum$Count<-table(cmzoo_df$SID)

ggdensity(data = cmzoo_df, x = "logKB", y = "..density..", fill = "SID",
          facet.by = "SID")
table(cmzoo_df$SID)
aggregate(cmzoo_df$KB ~ cmzoo_df$SID, FUN = "sum")

cm_cat2<-cm_cat
cm_cat2$Reference<-"rCheMyd1"
cm<-rbind(cmd_cat, cmzoo_cat, cm_cat2)
draftRef<-ggbarplot(cm, x = "Sample_ID", y = "sumMB", fill = "Category", orientation = "horiz",
          facet.by = "Reference")
ggsave(filename = "Plots/Inc_WGR/ROH_comparison_draft_DNAZoo.svg", draftRef,
       width = 29, height = 12, units = "cm")
