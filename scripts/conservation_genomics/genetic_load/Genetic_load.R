#########################################
### Genetic load summary after review ###
#########################################

library(ggpubr)
library(venn)

df<-read.csv(file="WGR/for_genome_paper/Genetic_load/snpEff_summary_forR.csv")

ggbarplot(df, x = "Sample", y = "Percentage", fill = "Spp", facet.by = "Variant_effect")

df$Percentage<-as.numeric(gsub("%","",df$Percentage))

high<-df[df$Variant_effect == "High",]
moderate<-df[df$Variant_effect == "Moderate",]
low<-df[df$Variant_effect == "Low",]
modifier<-df[df$Variant_effect == "Modifier",]

A<-ggbarplot(high, x = "Sample", y = "Percentage", fill = "Spp",
          title = "High Impact", ylab = "Percentage of variants") + rotate_x_text(45)
B<-ggbarplot(moderate, x = "Sample", y = "Percentage", fill = "Spp",
          title = "Moderate Impact", ylab = "Percentage of variants") + rotate_x_text(45)
C<-ggbarplot(low, x = "Sample", y = "Percentage", fill = "Spp",
          title = "Low Impact", ylab = "Percentage of variants") + rotate_x_text(45)
D<-ggbarplot(modifier, x = "Sample", y = "Percentage", fill = "Spp",
          title = "Modifier variants", ylab = "Percentage of variants") + rotate_x_text(45)
E<-ggarrange(A + rremove("xlab"),
             B + rremove("xlab"),
             C + rremove("xlab"),
             D + rremove("xlab"),
             ncol=4,nrow = 1, common.legend = T)
E
ggsave("Plots/Inc_WGR/genetic_load_FINAL.svg", E, width = 30, height = 10, units = "cm")

A<-ggboxplot(data = high, x = "Variant_effect", y = "Percentage", fill = "Spp",
             title = "High Impact", ylab = "Percentage of variants", xlab = "")
B<-ggboxplot(data = moderate, x = "Variant_effect", y = "Percentage", fill = "Spp",
             title = "Moderate Impact", ylab = "Percentage of variants", xlab = "")
C<-ggboxplot(data = low, x = "Variant_effect", y = "Percentage", fill = "Spp",
             title = "Low Impact", ylab = "Percentage of variants", xlab = "")
D<-ggboxplot(data = modifier, x = "Variant_effect", y = "Percentage", fill = "Spp",
             title = "Modifier variants", ylab = "Percentage of variants", xlab = "")

high_agg<-as.data.frame(merge(aggregate(high$Percentage ~ high$Spp, FUN = "mean"),
                        aggregate(high$Percentage ~ high$Spp, FUN = "sd"),
                        by = "high$Spp"))
colnames(high_agg)<-c("Species","Mean","SD")
high_agg$Upper<-high_agg$Mean+high_agg$SD
high_agg$Lower<-high_agg$Mean-high_agg$SD
A<-ggbarplot(high_agg, x = "Species", y = "Mean", fill = "Species",
          legend = "none",ylab = "Percentage of variants", xlab = "") + 
  geom_errorbar(aes(group = Species, ymax = Upper, ymin = Lower),
              position = position_dodge(width = 0.3), width = 0.05)

moderate_agg<-as.data.frame(merge(aggregate(moderate$Percentage ~ moderate$Spp, FUN = "mean"),
                              aggregate(moderate$Percentage ~ moderate$Spp, FUN = "sd"),
                              by = "moderate$Spp"))
colnames(moderate_agg)<-c("Species","Mean","SD")
moderate_agg$Upper<-moderate_agg$Mean+moderate_agg$SD
moderate_agg$Lower<-moderate_agg$Mean-moderate_agg$SD
B<-ggbarplot(moderate_agg, x = "Species", y = "Mean", fill = "Species",
          legend = "none",ylab = "Percentage of variants", xlab = "") + 
  geom_errorbar(aes(group = Species, ymax = Upper, ymin = Lower),
                position = position_dodge(width = 0.3), width = 0.05)

low_agg<-as.data.frame(merge(aggregate(low$Percentage ~ low$Spp, FUN = "mean"),
                              aggregate(low$Percentage ~ low$Spp, FUN = "sd"),
                              by = "low$Spp"))
colnames(low_agg)<-c("Species","Mean","SD")
low_agg$Upper<-low_agg$Mean+low_agg$SD
low_agg$Lower<-low_agg$Mean-low_agg$SD
C<-ggbarplot(low_agg, x = "Species", y = "Mean", fill = "Species",
          legend = "none",ylab = "Percentage of variants", xlab = "") + 
  geom_errorbar(aes(group = Species, ymax = Upper, ymin = Lower),
                position = position_dodge(width = 0.3), width = 0.05)

modifier_agg<-as.data.frame(merge(aggregate(modifier$Percentage ~ modifier$Spp, FUN = "mean"),
                              aggregate(modifier$Percentage ~ modifier$Spp, FUN = "sd"),
                              by = "modifier$Spp"))
colnames(modifier_agg)<-c("Species","Mean","SD")
modifier_agg$Upper<-modifier_agg$Mean+modifier_agg$SD
modifier_agg$Lower<-modifier_agg$Mean-modifier_agg$SD
D<-ggbarplot(modifier_agg, x = "Species", y = "Mean", fill = "Species",
          legend = "none",ylab = "Percentage of variants", xlab = "") + 
  geom_errorbar(aes(group = Species, ymax = Upper, ymin = Lower),
                position = position_dodge(width = 0.3), width = 0.05)

E<-ggarrange(A,B,C,D, ncol = 4, nrow = 1)

ggsave(filename = "Plots/Inc_WGR/Genetic_load_all_indv.svg", E,
       height = 12, width = 24, units = "cm")

t.test(high$Percentage ~ high$Spp)
t.test(moderate$Percentage ~ moderate$Spp)
t.test(low$Percentage ~ low$Spp)
t.test(modifier$Percentage ~ modifier$Spp)

## Mis:Silent
ratio<-as.data.frame(rbind(cbind(rbind(0.9816,0.984,0.9868,0.9811,0.9991),
                      rep("Dermochelys",5)),
              cbind(rbind(0.6959,0.7104,0.6972,0.6986,0.6985,0.6971),
                    rep("Chelonia",6))))
ratio$V1<-as.numeric(ratio$V1)
t.test(ratio$V1 ~ ratio$V2)

## By species:
dc<-df[df$Spp=="Dermochelys",]
cm<-df[df$Spp=="Chelonia",]

dch<-dc[dc$Variant_effect == "High",]
dcm<-dc[dc$Variant_effect == "Moderate",]
dcl<-dc[dc$Variant_effect == "Low",]
dcmo<-dc[dc$Variant_effect == "Modifier",]
D1<-ggbarplot(dch, x = "Sample", y = "Percentage", fill = "navy", title= 'High Impact')
D2<-ggbarplot(dcm, x = "Sample", y = "Percentage", fill = "navy", title= 'Moderate Impact')
D3<-ggbarplot(dcl, x = "Sample", y = "Percentage", fill = "navy", title= 'Low Impact')
D4<-ggbarplot(dcmo, x = "Sample", y = "Percentage", fill = "navy", title= 'Modifier')
ggarrange(D1, D2, D3, D4, ncol = 2, nrow = 2)

cmh<-cm[cm$Variant_effect == "High",]
cmm<-cm[cm$Variant_effect == "Moderate",]
cml<-cm[cm$Variant_effect == "Low",]
cmmo<-cm[cm$Variant_effect == "Modifier",]
C1<-ggbarplot(cmh, x = "Sample", y = "Percentage", fill = "dark red", title= 'High Impact') + rotate_x_text(45)
C2<-ggbarplot(cmm, x = "Sample", y = "Percentage", fill = "dark red", title= 'Moderate Impact')+ rotate_x_text(45)
C3<-ggbarplot(cml, x = "Sample", y = "Percentage", fill = "dark red", title= 'Low Impact')+ rotate_x_text(45)
C4<-ggbarplot(cmmo, x = "Sample", y = "Percentage", fill = "dark red", title= 'Modifier')+ rotate_x_text(45)
ggarrange(C1, C2, C3, C4, ncol = 2, nrow = 2)

### Gene lists:
# Dermochelys
dc_gene_list<-list.files(path = "WGR/for_genome_paper/Genetic_load/DerCor/", pattern = "txt")
dc_genes<-list()
for(q in 1:length(dc_gene_list)){
  df<-read.delim(file=paste0("WGR/for_genome_paper/Genetic_load/DerCor/",dc_gene_list[q]),
                 header = T, skip = 1, na.strings = c("","NA"))
  df$Sample<-gsub("_exons.genes.txt","", dc_gene_list[q])
  df$Spp<-"Dermochelys"
  df<-df[df$variants_impact_HIGH != 0,]
  df<-df[,1:8]
  dc_genes[[q]]<-df
}
dc_genes<-do.call(rbind, dc_genes)
dc_gene_names<-levels(factor(dc_genes$GeneId))
write.table(file = "WGR/for_genome_paper/Genetic_load/DerCor_gene_names.txt",
            dc_gene_names, quote = F, row.names = F, col.names = F)

# Chelonia
cm_gene_list<-list.files(path = "WGR/for_genome_paper/Genetic_load/CheMyd/", pattern = "txt")
cm_genes<-list()
for(q in 1:length(cm_gene_list)){
  df<-read.delim(file=paste0("WGR/for_genome_paper/Genetic_load/CheMyd/",cm_gene_list[q]),
                 header = T, skip = 1, na.strings = c("","NA"))
  df$Sample<-gsub("_exons.genes.txt","", cm_gene_list[q])
  df$Spp<-"Chelonia"
  df<-df[df$variants_impact_HIGH != 0,]
  cm_genes[[q]]<-df
}
cm_genes<-do.call(rbind, cm_genes)
cm_gene_names<-levels(factor(cm_genes$GeneId))
write.table(file = "WGR/for_genome_paper/Genetic_load/CheMyd_gene_names.txt",
            cm_gene_names, quote = F, row.names = F, col.names = F)


### Plot PANTHER
cm_BP<-read.delim(file="WGR/for_genome_paper/Genetic_load/CheMyd/pantherChart_CheMydHIGH.txt",
                  header = F)
cm_BP$Spp<-"Chelonia mydas"
dc_BP<-read.delim(file="WGR/for_genome_paper/Genetic_load/DerCor/pantherChart_DerCorHIGH.txt",
                  header = F)
dc_BP$Spp<-"Dermochelys coriacea"
BP<-as.data.frame(rbind(dc_BP, cm_BP))
BP$V4<-as.numeric(gsub("%","",BP$V4))
BP$V5<-as.numeric(gsub("%","",BP$V5))
A<-ggbarplot(BP, x = "V2", y = "V4", fill = "Spp", position = position_dodge(),
          orientation = "horiz", xlab = "", ylab = "Percent of genes")
ggsave(filename = "Plots/Inc_WGR/PANTHER_outs_high_impact_snpEff.svg", A,
       height = 21, width = 29, units = "cm")
