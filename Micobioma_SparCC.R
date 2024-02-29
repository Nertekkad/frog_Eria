#### Data of frog's micobiome under different temperatures ####

setwd("~/Frog's data")
# Data loading
library(readxl)
meta_hongos<-read.csv("Metadata_Hongos.csv")
hongos<-read_excel("Hongos.xlsx")

# Identify frog´s life stages
life_stage<-unique(meta_hongos$Life_stage)

# Separation of data by life stages type
Tadpole<-meta_hongos[which(meta_hongos$Life_stage == life_stage[1]),]
Metamorphic<-meta_hongos[which(meta_hongos$Life_stage == life_stage[2]),]
Sub_adult<-meta_hongos[which(meta_hongos$Life_stage == life_stage[3]),]

# Tadpole

# Treatment identification
Tad_T1<-Tadpole$SampleID[which(Tadpole$Treatment=="Treatment1")]
Tad_T2<-Tadpole$SampleID[which(Tadpole$Treatment=="Treatment2")]
Tad_Ctr<-Tadpole$SampleID[which(Tadpole$Treatment=="Control")]
# Separation of data by treatment types
HTad_T1<-as.data.frame(hongos[,which(colnames(hongos) %in% Tad_T1)])
HTad_T2<-as.data.frame(hongos[,which(colnames(hongos) %in% Tad_T2)])
HTad_Ctr<-as.data.frame(hongos[,which(colnames(hongos) %in% Tad_Ctr)])

# Metamorphic

# Treatment identification
Met_T1<-Metamorphic$SampleID[which(Metamorphic$Treatment=="Treatment1")]
Met_T2<-Metamorphic$SampleID[which(Metamorphic$Treatment=="Treatment2")]
Met_Ctr<-Metamorphic$SampleID[which(Metamorphic$Treatment=="Control")]
# Separation of data by treatment types
HMet_T1<-as.data.frame(hongos[,which(colnames(hongos) %in% Met_T1)])
HMet_T2<-as.data.frame(hongos[,which(colnames(hongos) %in% Met_T2)])
HMet_Ctr<-as.data.frame(hongos[,which(colnames(hongos) %in% Met_Ctr)])

# Sub-adult

# Treatment identification
Adl_T1<-Sub_adult$SampleID[which(Sub_adult$Treatment=="Treatment1")]
Adl_T2<-Sub_adult$SampleID[which(Sub_adult$Treatment=="Treatment2")]
Adl_Ctr<-Sub_adult$SampleID[which(Sub_adult$Treatment=="Control")]
# Separation of data by treatment types
HAdl_T1<-as.data.frame(hongos[,which(colnames(hongos) %in% Adl_T1)])
HAdl_T2<-as.data.frame(hongos[,which(colnames(hongos) %in% Adl_T2)])
HAdl_Ctr<-as.data.frame(hongos[,which(colnames(hongos) %in% Adl_Ctr)])

# Taxa table
t1<-which(colnames(hongos)=="Kingdom")
t2<-which(colnames(hongos)=="Species")
tax_fungi<-as.data.frame(hongos[,t1:t2])

##### Abundance tables' collapse #####

library(mlBioNets)

# Tadpole under treatment 1
FTad_T1<-T_collapse(F, T_table = tax_fungi, O_table = HTad_T1, names_level = "Genus")
FTad_T1<-FTad_T1[, -which(colnames(FTad_T1) == "unidentified")]
FTad_T1<-FTad_T1[, -which(is.na(colnames(FTad_T1)))]

# Tadpole under treatment 2
FTad_T2<-T_collapse(F, T_table = tax_fungi, O_table = HTad_T2, names_level = "Genus")
FTad_T2<-FTad_T2[, -which(colnames(FTad_T2) == "unidentified")]
FTad_T2<-FTad_T2[, -which(is.na(colnames(FTad_T2)))]

# Tadpole control
FTad_Ctr<-T_collapse(F, T_table = tax_fungi, O_table = HTad_Ctr, names_level = "Genus")
FTad_Ctr<-FTad_Ctr[, -which(colnames(FTad_Ctr) == "unidentified")]
FTad_Ctr<-FTad_Ctr[, -which(is.na(colnames(FTad_Ctr)))]

# Metamorphic under treatment 1
FMet_T1<-T_collapse(F, T_table = tax_fungi, O_table = HMet_T1, names_level = "Genus")
FMet_T1<-FMet_T1[, -which(colnames(FMet_T1) == "unidentified")]
FMet_T1<-FMet_T1[, -which(is.na(colnames(FMet_T1)))]

# Metamorphic under treatment 2
FMet_T2<-T_collapse(F, T_table = tax_fungi, O_table = HMet_T2, names_level = "Genus")
FMet_T2<-FMet_T2[, -which(colnames(FMet_T2) == "unidentified")]
FMet_T2<-FMet_T2[, -which(is.na(colnames(FMet_T2)))]

# Metamorphic control
FMet_Ctr<-T_collapse(F, T_table = tax_fungi, O_table = HMet_Ctr, names_level = "Genus")
FMet_Ctr<-FMet_Ctr[, -which(colnames(FMet_Ctr) == "unidentified")]
FMet_Ctr<-FMet_Ctr[, -which(is.na(colnames(FMet_Ctr)))]

# Sub-adult under treatment 1
FAdl_T1<-T_collapse(F, T_table = tax_fungi, O_table = HAdl_T1, names_level = "Genus")
FAdl_T1<-FAdl_T1[, -which(colnames(FAdl_T1) == "unidentified")]
FAdl_T1<-FAdl_T1[, -which(is.na(colnames(FAdl_T1)))]

# Sub-adult under treatment 2
FAdl_T2<-T_collapse(F, T_table = tax_fungi, O_table = HAdl_T2, names_level = "Genus")
FAdl_T2<-FAdl_T2[, -which(colnames(FAdl_T2) == "unidentified")]
FAdl_T2<-FAdl_T2[, -which(is.na(colnames(FAdl_T2)))]

# Sub-adult control
FAdl_Ctr<-T_collapse(F, T_table = tax_fungi, O_table = HAdl_Ctr, names_level = "Genus")
FAdl_Ctr<-FAdl_Ctr[, -which(colnames(FAdl_Ctr) == "unidentified")]
FAdl_Ctr<-FAdl_Ctr[, -which(is.na(colnames(FAdl_Ctr)))]

#### Network inference ####

##### ARACNe networks #####

# Tadpole
#FTad_T1Net<-net_inference(FTad_T1, "sparcc", 0.4)
#FTad_T2Net<-net_inference(FTad_T2, "sparcc", 0.4)
#FTad_CtrNet<-net_inference(FTad_Ctr, "sparcc", 0.4)
# Metamorphic
#FMet_T1Net<-net_inference(FMet_T1, "sparcc", 0.4)
#FMet_T2Net<-net_inference(FMet_T2, "sparcc", 0.4)
#FMet_CtrNet<-net_inference(FMet_Ctr, "sparcc", 0.4)
# Sub-adult
#FAdl_T1Net<-net_inference(FAdl_T1, "sparcc", 0.4)
#FAdl_T2Net<-net_inference(FAdl_T2, "sparcc", 0.4)
#FAdl_CtrNet<-net_inference(FAdl_Ctr, "sparcc", 0.4)

#ml_FTad_sp <- list(FTad_T2Net, FTad_T1Net, FTad_CtrNet)
#ml_FMet_sp <- list(FMet_T2Net, FMet_T1Net, FMet_CtrNet)
#ml_FAdl_sp <- list(FAdl_T2Net, FAdl_T1Net, FAdl_CtrNet)

#saveRDS(ml_FTad_sp, "~/frog_Eria/sparCC_Nets/ml_FTad_sp.RDS")
#saveRDS(ml_FMet_sp, "~/frog_Eria/sparCC_Nets/ml_FMet_sp.RDS")
#saveRDS(ml_FAdl_sp, "~/frog_Eria/sparCC_Nets/ml_FAdl_sp.RDS")

ml_FTad_sp<-readRDS("~/frog_Eria/sparCC_Nets/ml_FTad_sp.RDS")
ml_FMet_sp<-readRDS("~/frog_Eria/sparCC_Nets/ml_FMet_sp.RDS")
ml_FAdl_sp<-readRDS("~/frog_Eria/sparCC_Nets/ml_FAdl_sp.RDS")

#### Multilayer networks ####

library(viridis)
unq<-unique(tax_fungi[,"Phylum"])
unq<-unq[-c(which(unq == "unidentified"), which(is.na(unq)))]
colors <- sample(viridis(100), length(unq))

# Abundance tables list
abs_FTad<-list(FTad_T2, FTad_T1, FTad_Ctr) # Tadpole
abs_FMet<-list(FMet_T2, FMet_T1, FMet_Ctr) # Metamorphic
abs_FAdl<-list(FAdl_T2, FAdl_T1, FAdl_Ctr) # Adult

# Vertex colored by phylum
ml_FTad_sp <- v_colored_ml(ml_FTad_sp, tax_fungi, g_tax = "Phylum",
                           p_tax = "Genus", g_colors = colors)
ml_FMet_sp <- v_colored_ml(ml_FMet_sp, tax_fungi, g_tax = "Phylum",
                           p_tax = "Genus", g_colors = colors)
ml_FAdl_sp <- v_colored_ml(ml_FAdl_sp, tax_fungi, g_tax = "Phylum",
                           p_tax = "Genus", g_colors = colors)

plot(ml_FTad_sp[[1]], vertex.label.color="black",
     vertex.color = vertex.attributes(ml_FTad_sp[[1]])$color, vertex.label.cex=.5,
     vertex.label.dist=1,layout=layout_with_kk, vertex.size = 5,
     main = "Tadpole under treatment 1")
legend(x=-2.4, y=1, unq, title = "Mycobiome", pch=21, pt.bg=colors, pt.cex=1.3, cex=.8, bty="n", ncol=1)

library(muxViz)
library(mlBioNets)
# Tadpole
lay <- layoutMultiplex(ml_FTad_sp, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_FTad_sp, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_FTad, ml_FTad_sp, 10),
                 node.colors=node_color_mat(ml_FTad_sp, "phylo"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

# Metamorphic
lay <- layoutMultiplex(ml_FMet_sp, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_FMet_sp, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_FMet, ml_FMet_sp, 10),
                 node.colors=node_color_mat(ml_FMet_sp, "phylo"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

# Sub-adult
lay <- layoutMultiplex(ml_FAdl_sp, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_FAdl_sp, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_FAdl, ml_FAdl_sp, 10),
                 node.colors=node_color_mat(ml_FAdl_sp, "phylo"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

#### Centrality analysis ####

ml_FTad_sp <- ml_TaxGroup(ml_FTad_sp, tax_fungi, "Phylum", "Genus")
ml_FMet_sp <- ml_TaxGroup(ml_FMet_sp, tax_fungi, "Phylum", "Genus")
ml_FAdl_sp <- ml_TaxGroup(ml_FAdl_sp, tax_fungi, "Phylum", "Genus")

FTad_degree <- ctr_df(ml_FTad_sp, c("Treatment 2", "Treatment 1", "Control"))
FMet_degree <- ctr_df(ml_FMet_sp, c("Treatment 2", "Treatment 1", "Control"))
FAdl_degree <- ctr_df(ml_FAdl_sp, c("Treatment 2", "Treatment 1", "Control"))

FTad_phyl_degree <- phyl_ctr_df(FTad_degree, c("Treatment 2", "Treatment 1", "Control"),
                                n_layers = 3)
FMet_phyl_degree <- phyl_ctr_df(FMet_degree, c("Treatment 2", "Treatment 1", "Control"),
                                n_layers = 3)
FAdl_phyl_degree <- phyl_ctr_df(FAdl_degree, c("Treatment 2", "Treatment 1", "Control"),
                                n_layers = 3)

library(ggplot2)
ggplot(FTad_phyl_degree, aes(x = reorder(FTad_phyl_degree$Taxon, -FTad_phyl_degree$`Treatment 2`),
                             y = FTad_phyl_degree$`Treatment 2`, fill = FTad_phyl_degree$Taxon)) +
  geom_bar(stat = "identity") +
  labs(title = "Phyla importance by degree \n Treatment 2") +
  xlab("Phylum") + ylab("Degree") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(FTad_phyl_degree, aes(x = reorder(FTad_phyl_degree$Taxon, -FTad_phyl_degree$`Treatment 1`),
                             y = FTad_phyl_degree$`Treatment 1`, fill = FTad_phyl_degree$Taxon)) +
  geom_bar(stat = "identity") +
  labs(title = "Phyla importance by degree \n Treatment 1") +
  xlab("Phylum") + ylab("Degree") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(FTad_phyl_degree, aes(x = reorder(FTad_phyl_degree$Taxon, -FTad_phyl_degree$`Control`),
                             y = FTad_phyl_degree$`Control`, fill = FTad_phyl_degree$Taxon)) +
  geom_bar(stat = "identity") +
  labs(title = "Phyla importance by degree \n Control") +
  xlab("Phylum") + ylab("Degree") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##### Centrality log-fold change #####

###### Treatment 1 ######
log2fc <- -log2((FTad_phyl_degree$Control+1)/(FTad_phyl_degree$`Treatment 1`+1))
zscore <- (log2fc-mean(log2fc))/sd(log2fc)
df_degree <- data.frame(
  Phylum = FTad_phyl_degree$Taxon,
  log2fc = log2fc,
  z_score = zscore
)

library(ggpubr)

ggbarplot(df_degree, x = "Phylum", y = "z_score",
          fill = "Phylum",
          color = "white",
          palette = colors,
          sort.val = "desc",
          sort.by.groups = FALSE,
          x.text.angle = 90,
          ylab = "z_score",
          rotate = TRUE,
          ggtheme = theme_minimal()) +
  theme(legend.position = "none")

ggbarplot(df_degree, x = "Phylum", y = "log2fc",
          fill = "Phylum",
          color = "white",
          palette = colors,
          sort.val = "desc",
          sort.by.groups = FALSE,
          x.text.angle = 90,
          ylab = "log2fc",
          rotate = TRUE,
          ggtheme = theme_minimal()) +
  theme(legend.position = "none")
