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
FTad_T1Net<-net_inference(FTad_T1, "aracne")
FTad_T2Net<-net_inference(FTad_T2, "aracne")
FTad_CtrNet<-net_inference(FTad_Ctr, "aracne")
# Metamorphic
FMet_T1Net<-net_inference(FMet_T1, "aracne")
FMet_T2Net<-net_inference(FMet_T2, "aracne")
FMet_CtrNet<-net_inference(FMet_Ctr, "aracne")
# Sub-adult
FAdl_T1Net<-net_inference(FAdl_T1, "aracne")
FAdl_T2Net<-net_inference(FAdl_T2, "aracne")
FAdl_CtrNet<-net_inference(FAdl_Ctr, "aracne")

#### Multilayer networks ####

unq<-unique(tax_fungi[,"Phylum"])
unq<-unq[-c(which(unq == "unidentified"), which(is.na(unq)))]
colors <- sample(rainbow(100), length(unq))

# Life stages
ml_FTad<-list(FTad_T2Net, FTad_T1Net, FTad_CtrNet) # Tadpole
ml_FMet<-list(FMet_T2Net, FMet_T1Net, FMet_CtrNet) # Metamorphic
ml_FAdl<-list(FAdl_T2Net, FAdl_T1Net, FAdl_CtrNet) # Adult

# Abundance tables list
abs_FTad<-list(FTad_T2, FTad_T1, FTad_Ctr) # Tadpole
abs_FMet<-list(FMet_T2, FMet_T1, FMet_Ctr) # Metamorphic
abs_FAdl<-list(FAdl_T2, FAdl_T1, FAdl_Ctr) # Adult

# Vertex colored by phylum
ml_FTad <- v_colored_ml(ml_FTad, tax_fungi, g_tax = "Phylum",
                        p_tax = "Genus", g_colors = colors)
ml_FMet <- v_colored_ml(ml_FMet, tax_fungi, g_tax = "Phylum",
                        p_tax = "Genus", g_colors = colors)
ml_FAdl <- v_colored_ml(ml_FAdl, tax_fungi, g_tax = "Phylum",
                        p_tax = "Genus", g_colors = colors)

plot(ml_FTad[[1]], vertex.label.color="black",
     vertex.color = vertex.attributes(ml_FTad[[1]])$color, vertex.label.cex=.5,
     vertex.label.dist=1,layout=layout_with_kk, vertex.size = 5,
     main = "Tadpole under treatment 1")
legend(x=-2.4, y=1, unq, title = "Mycobiome", pch=21, pt.bg=colors, pt.cex=1.3, cex=.8, bty="n", ncol=1)


library(muxViz)
# Tadpole
lay <- layoutMultiplex(ml_FTad, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_FTad, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_FTad, ml_FTad, 10),
                 node.colors=node_color_mat(ml_FTad, "phylo"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

# Metamorphic
lay <- layoutMultiplex(ml_FMet, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_FMet, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_FMet, ml_FMet, 10),
                 node.colors=node_color_mat(ml_FMet, "phylo"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

# Sub-adult
lay <- layoutMultiplex(ml_FAdl, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_FAdl, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_FAdl, ml_FAdl, 10),
                 node.colors=node_color_mat(ml_FAdl, "phylo"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

#### Abundances distribtion ####

##### Rank-abundance analysis #####

# Tadpole
library(RADanalysis)
sample_classes <- c(rep(1, nrow(FTad_T2)),rep(2, nrow(FTad_T1)),
                    rep(3, nrow(FTad_Ctr)))
line_cols <- c("red3","darkorange1","green3")

FTad_mat<-rbind(FTad_T2, FTad_T1, FTad_Ctr)

# Normalized abundances
n_FTad_mat<-FTad_mat/rowSums(FTad_mat)

# Sort the abundances in a decreasing order
sorted_abs<-list()
for(i in 1:nrow(n_FTad_mat)){
  sorted_abs[[i]]<-sort(n_FTad_mat[i,], decreasing = T)
}
n_FTad_mat<-matrix(unlist(sorted_abs), nrow=length(sorted_abs),
                   ncol=length(sorted_abs[[1]]), byrow=TRUE)
#n_BTad_mat<-n_BTad_mat[,-which(colSums(n_BTad_mat)==0)]

# Plot the axis
plot(1e10,xlim = c(1,30),ylim = c(0,0.5),
     xlab = "Species rank",ylab = "Abundance",cex.lab = 1.5,axes = FALSE)
sfsmisc::eaxis(side = 1,at = c(1,5,10,15,20,25,30))
sfsmisc::eaxis(side = 2,at = c(0,0.1,0.2,0.3,0.4,0.5),las = 0)

# Rank-abundance colors per treatment
a <- representative_RAD(norm_rad = n_FTad_mat,sample_ids = which(sample_classes == 1),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[1],0.5),border = NA)
a <- representative_RAD(norm_rad = n_FTad_mat,sample_ids = which(sample_classes == 2),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[2],0.5),border = NA)
a <- representative_RAD(norm_rad = n_FTad_mat,sample_ids = which(sample_classes == 3),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[3],0.5),border = NA)
legend("topright",bty = "n",
       legend = c("Treatment 2","Treatment 1", "Control"),
       col = line_cols, lwd = 3, title = "Tadpole")




# Metamorphic
sample_classes <- c(rep(1, nrow(FMet_T2)),rep(2, nrow(FMet_T1)),
                    rep(3, nrow(FMet_Ctr)))
line_cols <- c("red3","darkorange1","green3")

FMet_mat<-rbind(FMet_T2, FMet_T1, FMet_Ctr)

# Normalized abundances
n_FMet_mat<-FMet_mat/rowSums(FMet_mat)

# Sort the abundances in a decreasing order
sorted_abs<-list()
for(i in 1:nrow(n_FMet_mat)){
  sorted_abs[[i]]<-sort(n_FMet_mat[i,], decreasing = T)
}
n_FMet_mat<-matrix(unlist(sorted_abs), nrow=length(sorted_abs),
                   ncol=length(sorted_abs[[1]]), byrow=TRUE)

# Plot the axis
plot(1e10,xlim = c(1,30),ylim = c(0,0.4),
     xlab = "Species rank",ylab = "Abundance",cex.lab = 1.5,axes = FALSE)
sfsmisc::eaxis(side = 1,at = c(1,5,10,15,20,25,30))
sfsmisc::eaxis(side = 2,at = c(0,0.1,0.2,0.3,0.4),las = 0)

# Rank-abundance colors per treatment
a <- representative_RAD(norm_rad = n_FMet_mat,sample_ids = which(sample_classes == 1),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[1],0.5),border = NA)
a <- representative_RAD(norm_rad = n_FMet_mat,sample_ids = which(sample_classes == 2),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[2],0.5),border = NA)
a <- representative_RAD(norm_rad = n_FMet_mat,sample_ids = which(sample_classes == 3),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[3],0.5),border = NA)
legend("topright",bty = "n",
       legend = c("Treatment 2","Treatment 1", "Control"),
       col = line_cols, lwd = 3, title = "Metamorphic")

# Adults
sample_classes <- c(rep(1, nrow(FAdl_T2)),rep(2, nrow(FAdl_T1)),
                    rep(3, nrow(FAdl_Ctr)))
line_cols <- c("red3","darkorange1","green3")

FAdl_mat<-rbind(FAdl_T1, FAdl_T2, FAdl_Ctr)

# Normalized abundances
n_FAdl_mat<-FAdl_mat/rowSums(FAdl_mat)

# Sort the abundances in a decreasing order
sorted_abs<-list()
for(i in 1:nrow(n_FAdl_mat)){
  sorted_abs[[i]]<-sort(n_FAdl_mat[i,], decreasing = T)
}
n_FAdl_mat<-matrix(unlist(sorted_abs), nrow=length(sorted_abs),
                   ncol=length(sorted_abs[[1]]), byrow=TRUE)

# Plot the axis
plot(1e10,xlim = c(1,30),ylim = c(0,0.4),
     xlab = "Species rank",ylab = "Abundance",cex.lab = 1.5,axes = FALSE)
sfsmisc::eaxis(side = 1,at = c(1,5,10,15,20,25,30))
sfsmisc::eaxis(side = 2,at = c(0,0.1,0.2,0.3,0.4),las = 0)

# Rank-abundance colors per treatment
a <- representative_RAD(norm_rad = n_FAdl_mat,sample_ids = which(sample_classes == 1),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[1],0.5),border = NA)
a <- representative_RAD(norm_rad = n_FAdl_mat,sample_ids = which(sample_classes == 2),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[2],0.5),border = NA)
a <- representative_RAD(norm_rad = n_FAdl_mat,sample_ids = which(sample_classes == 3),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[3],0.5),border = NA)
legend("topright",bty = "n",
       legend = c("Treatment 2","Treatment 1", "Control"),
       col = line_cols, lwd = 3, title = "Sub-adult")


# Tadpole
matplot(x = seq(1:ncol(n_FTad_mat)), y = t(n_FTad_mat),
        type = "l", xlab = "Species rank", ylab = "Abundance",
        main = "Log-log rank-abundance \n Tadpole", lwd = 1, lty=1,
        col=c(rep(line_cols[1], nrow(FTad_T2)),
              rep(line_cols[2], nrow(FTad_T1)),
              rep(line_cols[3], nrow(FTad_Ctr))),
        log = "xy")
# Metamorphic
matplot(x = seq(1:ncol(n_FMet_mat)), y = t(n_FMet_mat),
        type = "l", xlab = "Species rank", ylab = "Abundance",
        main = "Log-log rank-abundance \n Metamorphic", lwd = 1, lty=1,
        col=c(rep(line_cols[1], nrow(FMet_T2)),
              rep(line_cols[2], nrow(FMet_T1)),
              rep(line_cols[3], nrow(FMet_Ctr))),
        log = "xy")
# Sub-adult
matplot(x = seq(1:ncol(n_FAdl_mat)), y = t(n_FAdl_mat),
        type = "l", xlab = "Species rank", ylab = "Abundance",
        main = "Log-log rank-abundance \n Sub-adult", lwd = 1, lty=1,
        col=c(rep(line_cols[1], nrow(FAdl_T2)),
              rep(line_cols[2], nrow(FAdl_T1)),
              rep(line_cols[3], nrow(FAdl_Ctr))),
        log = "xy")

##### MSD analysis #####
# MSD analysis

# Tadpole
d <- dist(x = FTad_mat,method = "manhattan")
mds <- cmdscale(d = d,k = 5,eig = TRUE)
plot(mds$points,xlab = "First coordinate",ylab = "Second coordinate",pch = 19,
     cex =1,col = line_cols[sample_classes],
     main = "MDS analysis \n Tadpole")
sample_classes <- c(rep(1, nrow(FTad_T2)),rep(2, nrow(FTad_T1)),
                    rep(3, nrow(FTad_Ctr)))
a <- representative_point(input = mds$points,ids = which(sample_classes == 1),
                          col = scales::alpha(line_cols[1],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
a <- representative_point(input = mds$points,ids = which(sample_classes == 2),
                          col = scales::alpha(line_cols[2],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
a <- representative_point(input = mds$points,ids = which(sample_classes == 3),
                          col = scales::alpha(line_cols[3],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
legend("topright",bty = "n",legend = c("Treatment 2","Treatment 1","Control"),
       col = line_cols,pch = 19)

# Metamorphic
d <- dist(x = FMet_mat,method = "manhattan")
mds <- cmdscale(d = d,k = 5,eig = TRUE)
sample_classes <- c(rep(1, nrow(FMet_T2)),rep(2, nrow(FMet_T2)),
                    rep(3, nrow(FMet_Ctr)))
plot(mds$points,xlab = "First coordinate",ylab = "Second coordinate",pch = 19,
     cex =1,col = line_cols[sample_classes],
     main = "MDS analysis \n Metamorphic")
a <- representative_point(input = mds$points,ids = which(sample_classes == 1),
                          col = scales::alpha(line_cols[1],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
a <- representative_point(input = mds$points,ids = which(sample_classes == 2),
                          col = scales::alpha(line_cols[2],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
a <- representative_point(input = mds$points,ids = which(sample_classes == 3),
                          col = scales::alpha(line_cols[3],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
legend("bottomleft",bty = "n",legend = c("Treatment 2","Treatment 1","Control"),
       col = line_cols,pch = 19)

# Sub-adults
d <- dist(x = FAdl_mat,method = "manhattan")
mds <- cmdscale(d = d,k = 5,eig = TRUE)
sample_classes <- c(rep(1, nrow(FAdl_T2)),rep(2, nrow(FAdl_T1)),
                    rep(3, nrow(FAdl_Ctr)))
plot(mds$points,xlab = "First coordinate",ylab = "Second coordinate",pch = 19,
     cex =1,col = line_cols[sample_classes],
     main = "MDS analysis \n Sub-adults")
a <- representative_point(input = mds$points,ids = which(sample_classes == 1),
                          col = scales::alpha(line_cols[1],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
a <- representative_point(input = mds$points,ids = which(sample_classes == 2),
                          col = scales::alpha(line_cols[2],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
a <- representative_point(input = mds$points,ids = which(sample_classes == 3),
                          col = scales::alpha(line_cols[3],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
legend("bottomright",bty = "n",legend = c("Treatment 2","Treatment 1","Control"),
       col = line_cols,pch = 19)

#### Centrality analysis ####

ml_FTad <- ml_TaxGroup(ml_FTad, tax_fungi, "Phylum", "Genus")
ml_FMet <- ml_TaxGroup(ml_FMet, tax_fungi, "Phylum", "Genus")
ml_FAdl <- ml_TaxGroup(ml_FAdl, tax_fungi, "Phylum", "Genus")

FTad_degree <- ctr_df(ml_FTad, c("Treatment 2", "Treatment 1", "Control"))
FMet_degree <- ctr_df(ml_FMet, c("Treatment 2", "Treatment 1", "Control"))
FAdl_degree <- ctr_df(ml_FAdl, c("Treatment 2", "Treatment 1", "Control"))

FTad_phyl_degree <- phyl_ctr_df(FTad_degree, c("Treatment 2", "Treatment 1", "Control"),
                                n_layers = 3)
FMet_phyl_degree <- phyl_ctr_df(FMet_degree, c("Treatment 2", "Treatment 1", "Control"),
                                n_layers = 3)
FAdl_phyl_degree <- phyl_ctr_df(FAdl_degree, c("Treatment 2", "Treatment 1", "Control"),
                                n_layers = 3)

library(viridis)
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

###### Treatment 2 ######

log2fc <- -log2((FTad_phyl_degree$Control+1)/(FTad_phyl_degree$`Treatment 2`+1))
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




##### Degree distribution violin plot #####

T2<-degree(FTad_T2Net)[-which(degree(FTad_T2Net)==0)]
T1<-degree(FTad_T1Net)[-which(degree(FTad_T1Net)==0)]
Ctr<-degree(FTad_CtrNet)[-which(degree(FTad_CtrNet)==0)]

df_violinplot<-data.frame(
  Treatment = c(rep("Treatment 2", length(T2)), rep("Treatment 1", length(T1)),
                rep("Control", length(Ctr))),
  Data = c(T2, T1, Ctr)
)

library(ggstatsplot)
ggbetweenstats(
  data  = df_violinplot,
  x     = Treatment,
  y     = Data,
  title = "Tadpole"
)

T2<-degree(FMet_T2Net)[-which(degree(FMet_T2Net)==0)]
T1<-degree(FMet_T1Net)[-which(degree(FMet_T1Net)==0)]
Ctr<-degree(FMet_CtrNet)[-which(degree(FMet_CtrNet)==0)]

df_violinplot<-data.frame(
  Treatment = c(rep("Treatment 2", length(T2)), rep("Treatment 1", length(T1)),
                rep("Control", length(Ctr))),
  Data = c(T2, T1, Ctr)
)

library(ggstatsplot)
ggbetweenstats(
  data  = df_violinplot,
  x     = Treatment,
  y     = Data,
  title = "Metamorphic"
)

T2<-degree(FAdl_T2Net)[-which(degree(FAdl_T2Net)==0)]
T1<-degree(FAdl_T1Net)[-which(degree(FAdl_T1Net)==0)]
Ctr<-degree(FAdl_CtrNet)[-which(degree(FAdl_CtrNet)==0)]

df_violinplot<-data.frame(
  Treatment = c(rep("Treatment 2", length(T2)), rep("Treatment 1", length(T1)),
                rep("Control", length(Ctr))),
  Data = c(T2, T1, Ctr)
)

library(ggstatsplot)
ggbetweenstats(
  data  = df_violinplot,
  x     = Treatment,
  y     = Data,
  title = "Sub-adult"
)

##### Dysbiosis analysis #####

# Transform the data into a phyloseq object
tax_fungi<-as.matrix(hongos[,t1:t2])
otus_fungi<-as.matrix(hongos[,-(t1:t2)])
otus_fungi<-otus_fungi[,-1]
# OTUs IDs as row names and sample IDs as column names
ID_otus <- hongos$id
otus_fungi<-otus_fungi[,-1]
samples_fungi<-colnames(otus_fungi)
otus_fungi<-matrix(as.numeric(otus_fungi),
                   ncol = length(samples_fungi),
                   nrow = length(ID_otus))
colnames(otus_fungi)<-samples_fungi
rownames(otus_fungi)<-ID_otus
# OTUs IDs as row names in the taxa table
rownames(tax_fungi)<-ID_otus
# Phyloseq objects
library(phyloseq)
otus_fungi<-otu_table(otus_fungi, taxa_are_rows = TRUE)
tax_fungi<-tax_table(tax_fungi)
physeq = phyloseq(otus_fungi, tax_fungi)
# Sample metadata
rownames(meta_hongos)<-meta_hongos[,1]
meta_hongos<-meta_hongos[,-1]
sample_data<-sample_data(meta_hongos)
# Merge of phyloseq objects
physeq2<-merge_phyloseq(physeq, sample_data)

# Dysbiosis analysis
library(dysbiosisR)
# Bray-Curtis distance matrix
dist.mat <- phyloseq::distance(physeq2, "bray")
# Get reference samples
ref.samples <- sample_names(subset_samples(physeq2, 
                                           Treatment == "Control"))
# Community level variation analysis
dysbiosis_1 <- dysbiosisMedianCLV(physeq2,
                                  dist_mat = dist.mat,
                                  reference_samples = ref.samples)
# We sample the data set identifying as dysbiotic the data under the 90th percentile
dysbiosis_thres <- quantile(subset(dysbiosis_1, Treatment == "Treatment2")$score, 0.9)
normobiosis_thres <- quantile(subset(dysbiosis_1, Treatment == "Treatment2")$score, 0.1)

dysbiosis_1 <- dysbiosis_1 |> 
  mutate(isDysbiostic = ifelse(score >= dysbiosis_thres, TRUE, FALSE))

# Dysbiosis plot measures according to CLV method
p1 <- plotDysbiosis(df=dysbiosis_1,
                    xvar="Treatment",
                    yvar="score",
                    colors=c(Treatment1="orange", Treatment2="red",
                             Control="green"),
                    show_points = FALSE) +
  labs(x="", y="Dysbiosis Score") +
  theme_bw(base_size = 14)
p1

# Dysbiosis plot measures according to euclidean method
dysbiosis_2 <- euclideanDistCentroids(physeq2,
                                      dist_mat = dist.mat,
                                      use_squared = TRUE,
                                      group_col = "Treatment",
                                      control_label = "Control",
                                      case_label = "Treatment2")

p2 <- plotDysbiosis(df=dysbiosis_2,
                    xvar="Treatment",
                    yvar="CentroidDist_score",
                    colors=c(Treatment1="orange", Treatment2="red",
                             Control="green"),
                    show_points = FALSE) +
  labs(x="", y="Dysbiosis Score (Centroid)") +
  theme_bw(base_size = 14)
p2


# Dysbiosis plot measures according to the Combined alpha-beta diversity based score
dysbiosis_3 <- combinedShannonJSD(physeq2,
                                  reference_samples = ref.samples)


p3 <- plotDysbiosis(df=dysbiosis_3,
                    xvar="Treatment",
                    yvar="ShannonJSDScore",
                    colors=c(Treatment1="orange", Treatment2="red",
                             Control="green"),
                    show_points = FALSE)+
  labs(x="", y="Shannon-JSD\nDysbiosis Score") +
  theme_bw(base_size = 14)
p3


cloud.results <- cloudStatistic(physeq2,
                                dist_mat = dist.mat,
                                reference_samples = ref.samples,
                                ndim=-1,
                                k_num=5)

p4 <- plotDysbiosis(df=cloud.results,
                    xvar="Treatment",
                    yvar="log2Stats",
                    colors=c(Treatment1="orange", Treatment2="red",
                             Control="green"),
                    show_points = FALSE) +
  labs(x="", y="Dysbiosis CLOUD Score") +
  theme_bw(base_size = 14)
p4

