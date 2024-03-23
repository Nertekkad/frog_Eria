#### Data of frog's bacteriome under different temperatures ####

setwd("~/Frog's data")
# Data loading
library(readxl)
meta_bacterias<-read.csv("Metadata_Bacterias.csv")
bacterias<-read_excel("Bacterias.xlsx")

# Identify frog´s life stages
life_stage<-unique(meta_bacterias$Life_stage)

##### Separation of data by life stages type #####
Tadpole<-meta_bacterias[which(meta_bacterias$Life_stage == life_stage[1]),]
Metamorphic<-meta_bacterias[which(meta_bacterias$Life_stage == life_stage[2]),]
Sub_adult<-meta_bacterias[which(meta_bacterias$Life_stage == life_stage[3]),]

# Tadpole

# Treatment identification
Tad_T1<-Tadpole$sample.id[which(Tadpole$Treatment=="Treatment1")]
Tad_T2<-Tadpole$sample.id[which(Tadpole$Treatment=="Treatment2")]
Tad_Ctr<-Tadpole$sample.id[which(Tadpole$Treatment=="Control")]
# Separation of data by treatment types
BaTad_T1<-as.data.frame(bacterias[,which(colnames(bacterias) %in% Tad_T1)])
BaTad_T2<-as.data.frame(bacterias[,which(colnames(bacterias) %in% Tad_T2)])
BaTad_Ctr<-as.data.frame(bacterias[,which(colnames(bacterias) %in% Tad_Ctr)])

# Metamorphic

# Treatment identification
Met_T1<-Metamorphic$sample.id[which(Metamorphic$Treatment=="Treatment1")]
Met_T2<-Metamorphic$sample.id[which(Metamorphic$Treatment=="Treatment2")]
Met_Ctr<-Metamorphic$sample.id[which(Metamorphic$Treatment=="Control")]
# Separation of data by treatment types
BaMet_T1<-as.data.frame(bacterias[,which(colnames(bacterias) %in% Met_T1)])
BaMet_T2<-as.data.frame(bacterias[,which(colnames(bacterias) %in% Met_T2)])
BaMet_Ctr<-as.data.frame(bacterias[,which(colnames(bacterias) %in% Met_Ctr)])

# Sub-adult

# Treatment identification
Adl_T1<-Sub_adult$sample.id[which(Sub_adult$Treatment=="Treatment1")]
Adl_T2<-Sub_adult$sample.id[which(Sub_adult$Treatment=="Treatment2")]
Adl_Ctr<-Sub_adult$sample.id[which(Sub_adult$Treatment=="Control")]
# Separation of data by treatment types
BaAdl_T1<-as.data.frame(bacterias[,which(colnames(bacterias) %in% Adl_T1)])
BaAdl_T2<-as.data.frame(bacterias[,which(colnames(bacterias) %in% Adl_T2)])
BaAdl_Ctr<-as.data.frame(bacterias[,which(colnames(bacterias) %in% Adl_Ctr)])

# Taxa table
t1<-which(colnames(bacterias)=="Kingdom")
t2<-which(colnames(bacterias)=="Species")
tax_bacter<-as.data.frame(bacterias[,t1:t2])


##### Abundance tables' collapse #####

library(mlBioNets)

# Tadpole under treatment 1
BTad_T1<-T_collapse(F, T_table = tax_bacter, O_table = BaTad_T1, names_level = "Genus")
BTad_T1<-BTad_T1[, -which(colnames(BTad_T1) == "uncultured")]
BTad_T1<-BTad_T1[, -which(is.na(colnames(BTad_T1)))]

# Tadpole under treatment 2
BTad_T2<-T_collapse(F, T_table = tax_bacter, O_table = BaTad_T2, names_level = "Genus")
BTad_T2<-BTad_T2[, -which(colnames(BTad_T2) == "uncultured")]
BTad_T2<-BTad_T2[, -which(is.na(colnames(BTad_T2)))]

# Tadpole control
BTad_Ctr<-T_collapse(F, T_table = tax_bacter, O_table = BaTad_Ctr, names_level = "Genus")
BTad_Ctr<-BTad_Ctr[, -which(colnames(BTad_Ctr) == "uncultured")]
BTad_Ctr<-BTad_Ctr[, -which(is.na(colnames(BTad_Ctr)))]

# Metamorphic under treatment 1
BMet_T1<-T_collapse(F, T_table = tax_bacter, O_table = BaMet_T1, names_level = "Genus")
BMet_T1<-BMet_T1[, -which(colnames(BMet_T1) == "uncultured")]
BMet_T1<-BMet_T1[, -which(is.na(colnames(BMet_T1)))]

# Metamorphic under treatment 2
BMet_T2<-T_collapse(F, T_table = tax_bacter, O_table = BaMet_T2, names_level = "Genus")
BMet_T2<-BMet_T2[, -which(colnames(BMet_T2) == "uncultured")]
BMet_T2<-BMet_T2[, -which(is.na(colnames(BMet_T2)))]

# Metamorphic control
BMet_Ctr<-T_collapse(F, T_table = tax_bacter, O_table = BaMet_Ctr, names_level = "Genus")
BMet_Ctr<-BMet_Ctr[, -which(colnames(BMet_Ctr) == "uncultured")]
BMet_Ctr<-BMet_Ctr[, -which(is.na(colnames(BMet_Ctr)))]

# Sub-adult under treatment 1
BAdl_T1<-T_collapse(F, T_table = tax_bacter, O_table = BaAdl_T1, names_level = "Genus")
BAdl_T1<-BAdl_T1[, -which(colnames(BAdl_T1) == "uncultured")]
BAdl_T1<-BAdl_T1[, -which(is.na(colnames(BAdl_T1)))]

# Sub-adult under treatment 2
BAdl_T2<-T_collapse(F, T_table = tax_bacter, O_table = BaAdl_T2, names_level = "Genus")
BAdl_T2<-BAdl_T2[, -which(colnames(BAdl_T2) == "uncultured")]
BAdl_T2<-BAdl_T2[, -which(is.na(colnames(BAdl_T2)))]

# Sub-adult control
BAdl_Ctr<-T_collapse(F, T_table = tax_bacter, O_table = BaAdl_T2, names_level = "Genus")
BAdl_Ctr<-BAdl_Ctr[, -which(colnames(BAdl_Ctr) == "uncultured")]
BAdl_Ctr<-BAdl_Ctr[, -which(is.na(colnames(BAdl_Ctr)))]


#### Network inference ####

##### ARACNe networks #####

# Tadpole
#BTad_T1Net<-net_inference(BTad_T1, "sparcc", 0.4)
#BTad_T2Net<-net_inference(BTad_T2, "sparcc", 0.4)
#BTad_CtrNet<-net_inference(BTad_Ctr, "sparcc", 0.4)
# Metamorphic
#BMet_T1Net<-net_inference(BMet_T1, "sparcc", 0.4)
#BMet_T2Net<-net_inference(BMet_T2, "sparcc", 0.4)
#BMet_CtrNet<-net_inference(BMet_Ctr, "sparcc", 0.4)
# Sub-adult
#BAdl_T1Net<-net_inference(BAdl_T1, "sparcc", 0.4)
#BAdl_T2Net<-net_inference(BAdl_T2, "sparcc", 0.4)
#BAdl_CtrNet<-net_inference(BAdl_Ctr, "sparcc", 0.4)

ml_BTad_sp<-readRDS("~/frog_Eria/sparCC_Nets/ml_BTad2.RDS")
ml_BMet_sp<-readRDS("~/frog_Eria/sparCC_Nets/ml_BMet2.RDS")
ml_BAdl_sp<-readRDS("~/frog_Eria/sparCC_Nets/ml_BAdl2.RDS")

#### Multilayer networks ####

unq<-unique(tax_bacter[,"Phylum"])
unq<-unq[-c(which(unq == "uncultured"), which(is.na(unq)))]
colors <- sample(viridis(100), length(unq))

# Abundance tables list
abs_BTad<-list(BTad_T2, BTad_T1, BTad_Ctr) # Tadpole
abs_BMet<-list(BMet_T2, BMet_T1, BMet_Ctr) # Metamorphic
abs_BAdl<-list(BAdl_T2, BAdl_T1, BAdl_Ctr) # Adult

ml_BTad_sp <- v_colored_ml(ml_BTad_sp, tax_bacter, g_tax = "Phylum",
                         p_tax = "Genus", g_colors = colors)
ml_BMet_sp <- v_colored_ml(ml_BMet_sp, tax_bacter, g_tax = "Phylum",
                         p_tax = "Genus", g_colors = colors)
ml_BAdl_sp <- v_colored_ml(ml_BAdl_sp, tax_bacter, g_tax = "Phylum",
                         p_tax = "Genus", g_colors = colors)


library(igraph)
plot(ml_BTad_sp[[1]], vertex.label.color="black",
     vertex.color = vertex.attributes(ml_BTad_sp[[1]])$color, vertex.label.cex=.5,
     vertex.label.dist=1,layout=layout_with_kk, vertex.size = 5,
     main = "Tadpole under treatment 1")
legend(x=-2.4, y=1.4, unq, title = "Bacteriome", pch=21, pt.bg=colors, pt.cex=1.3, cex=.8, bty="n", ncol=1)

library(muxViz)
# Tadpole
lay <- layoutMultiplex(ml_BTad_sp, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_BTad_sp, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_BTad, ml_BTad_sp, 2),
                 node.colors=node_color_mat(ml_BTad_sp, "phylo"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)
# Metamorphic
lay <- layoutMultiplex(ml_BMet_sp, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_BMet_sp, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_BMet, ml_BMet_sp, 2),
                 node.colors=node_color_mat(ml_BMet_sp, "phylo"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)
# Sub-adult
lay <- layoutMultiplex(ml_BAdl_sp, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_BAdl_sp, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_BAdl, ml_BAdl_sp, 2),
                 node.colors=node_color_mat(ml_BAdl_sp, "phylo"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)


#### Centrality analysis ####

ml_BTad_sp <- ml_TaxGroup(ml_BTad_sp, tax_bacter, "Phylum", "Genus")
ml_BMet_sp <- ml_TaxGroup(ml_BMet_sp, tax_bacter, "Phylum", "Genus")
ml_BAdl_sp <- ml_TaxGroup(ml_BAdl_sp, tax_bacter, "Phylum", "Genus")

BTad_degree <- ctr_df(ml_BTad_sp, c("Treatment 2", "Treatment 1", "Control"))
BMet_degree <- ctr_df(ml_BMet_sp, c("Treatment 2", "Treatment 1", "Control"))
BAdl_degree <- ctr_df(ml_BAdl_sp, c("Treatment 2", "Treatment 1", "Control"))

BTad_phyl_degree <- phyl_ctr_df(BTad_degree, c("Treatment 2", "Treatment 1", "Control"),
                                n_layers = 3)
BMet_phyl_degree <- phyl_ctr_df(BMet_degree, c("Treatment 2", "Treatment 1", "Control"),
                                n_layers = 3)
BAdl_phyl_degree <- phyl_ctr_df(BAdl_degree, c("Treatment 2", "Treatment 1", "Control"),
                                n_layers = 3)

library(gridExtra)
grid.arrange(
  phyl_ctr_plot(BTad_phyl_degree, "Treatment 2", "Phyla importance by degree \n Treatment 2"),
  phyl_ctr_plot(BTad_phyl_degree, "Treatment 1", "Phyla importance by degree \n Treatment 1"),
  phyl_ctr_plot(BTad_phyl_degree, "Control", "Phyla importance by degree \n Control")
)

grid.arrange(
  phyl_ctr_plot(BMet_phyl_degree, "Treatment 2", "Phyla importance by degree \n Treatment 2"),
  phyl_ctr_plot(BMet_phyl_degree, "Treatment 1", "Phyla importance by degree \n Treatment 1"),
  phyl_ctr_plot(BMet_phyl_degree, "Control", "Phyla importance by degree \n Control")
)

grid.arrange(
  phyl_ctr_plot(BAdl_phyl_degree, "Treatment 2", "Phyla importance by degree \n Treatment 2"),
  phyl_ctr_plot(BAdl_phyl_degree, "Treatment 1", "Phyla importance by degree \n Treatment 1"),
  phyl_ctr_plot(BAdl_phyl_degree, "Control", "Phyla importance by degree \n Control")
)


##### Centrality log-fold change #####


##### Centrality log-fold change #####

library(gridExtra)
grid.arrange(log2fc(BTad_phyl_degree, "Control", "Treatment 1"),
             log2fc(BTad_phyl_degree, "Control", "Treatment 2"),
             log2fc(BTad_phyl_degree, "Treatment 1", "Treatment 2"),
             ncol = 3)


grid.arrange(log2fc(BMet_phyl_degree, "Control", "Treatment 1"),
             log2fc(BMet_phyl_degree, "Control", "Treatment 2"),
             log2fc(BMet_phyl_degree, "Treatment 1", "Treatment 2"),
             ncol = 3)

grid.arrange(log2fc(BAdl_phyl_degree, "Control", "Treatment 1"),
             log2fc(BAdl_phyl_degree, "Control", "Treatment 2"),
             log2fc(BAdl_phyl_degree, "Treatment 1", "Treatment 2"),
             ncol = 3)

##### Multiplex centrality #####

###### Degree ######

ml_BTad_degree <- ctr_ml(ml_BTad_sp, "degree")
ml_BMet_degree <- ctr_ml(ml_BMet_sp, "degree")
ml_BAdl_degree <- ctr_ml(ml_BAdl_sp, "degree")

# Tadpole
lay <- layoutMultiplex(ml_BTad_degree, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_BTad_degree, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_BTad, ml_BTad_degree, 2),
                 node.colors=node_color_mat(ml_BTad_degree, "centrality"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

# Metamorphic
lay <- layoutMultiplex(ml_BMet_degree, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_BMet_degree, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_BMet, ml_BMet_degree, 2),
                 node.colors=node_color_mat(ml_BMet_degree, "centrality"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

# Sub-adult
lay <- layoutMultiplex(ml_BAdl_degree, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_BAdl_degree, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_BAdl, ml_BAdl_degree, 2),
                 node.colors=node_color_mat(ml_BAdl_degree, "centrality"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)


#### Abundance vs degree ####

layer_names<-c("Control", "Treatment 1", "Treatment 2")
layer_colors<-c("green", "orange", "red")

cor_degree_abs(ml_BTad_sp, abs_BTad, layer_names, layer_colors, "Tadpole")
cor_degree_abs(ml_BMet_sp, abs_BMet, layer_names, layer_colors, "Metamorphic")
cor_degree_abs(ml_BAdl_sp, abs_BAdl, layer_names, layer_colors, "Sub-adult")


#### Topological properties ####

treatments <- c(20, 25, 29)
frog_ml_properties<-rbind(ml_properties(ml_BTad_sp, treatments),
                          ml_properties(ml_BMet_sp, treatments),
                          ml_properties(ml_BAdl_sp, treatments))

Stage = c(rep("Tadpole", length(ml_BTad_sp)),
          rep("Metamorphic",length(ml_BTad_sp)),
          rep("Sub-adult", length(ml_BTad_sp)))

frog_ml_properties <- cbind(frog_ml_properties, Stage)

require(ggplot2)

p1<-ggplot(frog_ml_properties, aes(x = Treatments, y = Mean_degree, color = Stage)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  
  labs(title = "",
       x = "Temperature",
       y = "Mean degree") +
  scale_color_manual(values = c("darkgreen", "purple4", "brown4")) +
  theme_minimal()

p2<-ggplot(frog_ml_properties, aes(x = Treatments, y = sd_degree, color = Stage)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  
  labs(title = "",
       x = "Temperature",
       y = "SD of degree") +
  scale_color_manual(values = c("darkgreen", "purple4", "brown4")) +
  theme_minimal()

p3<-ggplot(frog_ml_properties, aes(x = Treatments, y = Clusterization, color = Stage)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  
  labs(title = "",
       x = "Temperature",
       y = "Transitivity") +
  scale_color_manual(values = c("darkgreen", "purple4", "brown4")) +
  theme_minimal()

p4<-ggplot(frog_ml_properties, aes(x = Treatments, y = Edge_density, color = Stage)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  
  labs(title = "",
       x = "Temperature",
       y = "Edge density") +
  scale_color_manual(values = c("darkgreen", "purple4", "brown4")) +
  theme_minimal()

p5<-ggplot(frog_ml_properties, aes(x = Treatments, y = Connected_nodes, color = Stage)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  
  labs(title = "",
       x = "Temperature",
       y = "Proportion of linked nodes") +
  scale_color_manual(values = c("darkgreen", "purple4", "brown4")) +
  theme_minimal()

p6<-ggplot(frog_ml_properties, aes(x = Treatments, y = Modularity, color = Stage)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  
  labs(title = "",
       x = "Temperature",
       y = "Modularity") +
  scale_color_manual(values = c("darkgreen", "purple4", "brown4")) +
  theme_minimal()

grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3)

# Table

library(kableExtra)
kable(frog_ml_properties) %>% kable_styling()

# Slopes and R2

kable(slope_R2(ml_properties(ml_BTad_sp, treatments))) %>% kable_styling()
kable(slope_R2(ml_properties(ml_BMet_sp, treatments))) %>% kable_styling()
kable(slope_R2(ml_properties(ml_BAdl_sp, treatments))) %>% kable_styling()
