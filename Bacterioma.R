#### Data of frog's bacteriome under different temperatures ####

# Data of frog's microbiota under different temperatures

setwd("~/Frog's data")
# Data loading
library(readxl)
meta_bacterias<-read.csv("Metadata_Bacterias.csv")
bacterias<-read_excel("Bacterias.xlsx")

# Identify frog´s life stages
life_stage<-unique(meta_bacterias$Life_stage)

# Separation of data by life stages type
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

# Abundance tables' collapse

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
BAdl_Ctr<-T_collapse(F, T_table = tax_bacter, O_table = BaAdl_Ctr, names_level = "Genus")
BAdl_Ctr<-BAdl_Ctr[, -which(colnames(BAdl_Ctr) == "uncultured")]
BAdl_Ctr<-BAdl_Ctr[, -which(is.na(colnames(BAdl_Ctr)))]


#### Network inference ####

##### ARACNe networks #####

# Tadpole
BTad_T1Net<-net_inference(BTad_T1, "aracne")
BTad_T2Net<-net_inference(BTad_T2, "aracne")
BTad_CtrNet<-net_inference(BTad_Ctr, "aracne")
# Metamorphic
BMet_T1Net<-net_inference(BMet_T1, "aracne")
BMet_T2Net<-net_inference(BMet_T2, "aracne")
BMet_CtrNet<-net_inference(BMet_Ctr, "aracne")
# Sub-adult
BAdl_T1Net<-net_inference(BAdl_T1, "aracne")
BAdl_T2Net<-net_inference(BAdl_T2, "aracne")
BAdl_CtrNet<-net_inference(BAdl_Ctr, "aracne")


#### Multilayer networks ####

library(viridis)
unq<-unique(tax_bacter[,"Phylum"])
unq<-unq[-c(which(unq == "uncultured"), which(is.na(unq)))]
colors <- sample(viridis(100), length(unq))

# Life stages
ml_BTad<-list(BTad_T2Net, BTad_T1Net, BTad_CtrNet) # Tadpole
ml_BMet<-list(BMet_T2Net, BMet_T1Net, BMet_CtrNet) # Metamorphic
ml_BAdl<-list(BAdl_T2Net, BAdl_T1Net, BAdl_CtrNet) # Adult

# Abundance tables list
abs_BTad<-list(BTad_T2, BTad_T1, BTad_Ctr) # Tadpole
abs_BMet<-list(BMet_T2, BMet_T1, BMet_Ctr) # Metamorphic
abs_BAdl<-list(BAdl_T2, BAdl_T1, BAdl_Ctr) # Adult

# Vertex colored by phylum
ml_BTad <- v_colored_ml(ml_BTad, tax_bacter, g_tax = "Phylum",
                        p_tax = "Genus", g_colors = colors)
ml_BMet <- v_colored_ml(ml_BMet, tax_bacter, g_tax = "Phylum",
                        p_tax = "Genus", g_colors = colors)
ml_BAdl <- v_colored_ml(ml_BAdl, tax_bacter, g_tax = "Phylum",
                        p_tax = "Genus", g_colors = colors)

plot(ml_BTad[[1]], vertex.label.color="black",
     vertex.color = vertex.attributes(BTad_T1Net)$color, vertex.label.cex=.5,
     vertex.label.dist=1,layout=layout_with_kk, vertex.size = 5,
     main = "Tadpole under treatment 1")
legend(x=-2.4, y=1.4, unq, title = "Bacteriome", pch=21, pt.bg=colors, pt.cex=1.3, cex=.8, bty="n", ncol=1)


library(muxViz)
# Tadpole
lay <- layoutMultiplex(ml_BTad, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_BTad, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_BTad, ml_BTad, 2),
                 node.colors=node_color_mat(ml_BTad, "phylo"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)
# Metamorphic
lay <- layoutMultiplex(ml_BMet, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_BMet, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_BMet, ml_BMet, 2),
                 node.colors=node_color_mat(ml_BMet, "phylo"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)
# Sub-adult
lay <- layoutMultiplex(ml_BAdl, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_BAdl, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_BAdl, ml_BAdl, 2),
                 node.colors=node_color_mat(ml_BAdl, "phylo"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

#### Abundances distribtion ####

##### Rank-abundance analysis #####

# Tadpole
library(RADanalysis)
sample_classes <- c(rep(1, nrow(BTad_T2)),rep(2, nrow(BTad_T1)),
                    rep(3, nrow(BTad_Ctr)))
line_cols <- c("red3","darkorange1","green3")

BTad_mat<-rbind(BTad_T2, BTad_T1, BTad_Ctr)

# Normalized abundances
n_BTad_mat<-BTad_mat/rowSums(BTad_mat)

# Sort the abundances in a decreasing order
sorted_abs<-list()
for(i in 1:nrow(n_BTad_mat)){
  sorted_abs[[i]]<-sort(n_BTad_mat[i,], decreasing = T)
}
n_BTad_mat<-matrix(unlist(sorted_abs), nrow=length(sorted_abs),
                   ncol=length(sorted_abs[[1]]), byrow=TRUE)
#n_BTad_mat<-n_BTad_mat[,-which(colSums(n_BTad_mat)==0)]

# Plot the axis
plot(1e10,xlim = c(1,30),ylim = c(0,0.4),
     xlab = "Species rank",ylab = "Abundance",cex.lab = 1.5,axes = FALSE)
sfsmisc::eaxis(side = 1,at = c(1,5,10,15,20,25,30))
sfsmisc::eaxis(side = 2,at = c(0,0.1,0.2,0.3,0.4),las = 0)

# Rank-abundance colors per treatment
a <- representative_RAD(norm_rad = n_BTad_mat,sample_ids = which(sample_classes == 1),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[1],0.5),border = NA)
a <- representative_RAD(norm_rad = n_BTad_mat,sample_ids = which(sample_classes == 2),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[2],0.5),border = NA)
a <- representative_RAD(norm_rad = n_BTad_mat,sample_ids = which(sample_classes == 3),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[3],0.5),border = NA)
legend("topright",bty = "n",
       legend = c("Treatment 2","Treatment 1", "Control"),
       col = line_cols, lwd = 3, title = "Tadpole")




# Metamorphic
sample_classes <- c(rep(1, nrow(BMet_T2)),rep(2, nrow(BMet_T1)),
                    rep(3, nrow(BMet_Ctr)))
line_cols <- c("red3","darkorange1","green3")

BMet_mat<-rbind(BMet_T2, BMet_T1, BMet_Ctr)

# Normalized abundances
n_BMet_mat<-BMet_mat/rowSums(BMet_mat)

# Sort the abundances in a decreasing order
sorted_abs<-list()
for(i in 1:nrow(n_BMet_mat)){
  sorted_abs[[i]]<-sort(n_BMet_mat[i,], decreasing = T)
}
n_BMet_mat<-matrix(unlist(sorted_abs), nrow=length(sorted_abs),
                   ncol=length(sorted_abs[[1]]), byrow=TRUE)

# Plot the axis
plot(1e10,xlim = c(1,30),ylim = c(0,0.3),
     xlab = "Species rank",ylab = "Abundance",cex.lab = 1.5,axes = FALSE)
sfsmisc::eaxis(side = 1,at = c(1,5,10,15,20,25,30))
sfsmisc::eaxis(side = 2,at = c(0,0.1,0.2,0.3),las = 0)

# Rank-abundance colors per treatment
a <- representative_RAD(norm_rad = n_BMet_mat,sample_ids = which(sample_classes == 1),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[1],0.5),border = NA)
a <- representative_RAD(norm_rad = n_BMet_mat,sample_ids = which(sample_classes == 2),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[2],0.5),border = NA)
a <- representative_RAD(norm_rad = n_BMet_mat,sample_ids = which(sample_classes == 3),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[3],0.5),border = NA)
legend("topright",bty = "n",
       legend = c("Treatment 2","Treatment 1", "Control"),
       col = line_cols, lwd = 3, title = "Metamorphic")

# Adults
sample_classes <- c(rep(1, nrow(BAdl_T2)),rep(2, nrow(BAdl_T1)),
                    rep(3, nrow(BAdl_Ctr)))
line_cols <- c("red3","darkorange1","green3")

BAdl_mat<-rbind(BAdl_T1, BAdl_T2, BAdl_Ctr)

# Normalized abundances
n_BAdl_mat<-BAdl_mat/rowSums(BAdl_mat)

# Sort the abundances in a decreasing order
sorted_abs<-list()
for(i in 1:nrow(n_BAdl_mat)){
  sorted_abs[[i]]<-sort(n_BAdl_mat[i,], decreasing = T)
}
n_BAdl_mat<-matrix(unlist(sorted_abs), nrow=length(sorted_abs),
                   ncol=length(sorted_abs[[1]]), byrow=TRUE)

# Plot the axis
plot(1e10,xlim = c(1,30),ylim = c(0,0.6),
     xlab = "Species rank",ylab = "Abundance",cex.lab = 1.5,axes = FALSE)
sfsmisc::eaxis(side = 1,at = c(1,5,10,15,20,25,30))
sfsmisc::eaxis(side = 2,at = c(0,0.1,0.2,0.3,0.4,0.5,0.6),las = 0)

# Rank-abundance colors per treatment
a <- representative_RAD(norm_rad = n_BAdl_mat,sample_ids = which(sample_classes == 1),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[1],0.5),border = NA)
a <- representative_RAD(norm_rad = n_BAdl_mat,sample_ids = which(sample_classes == 2),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[2],0.5),border = NA)
a <- representative_RAD(norm_rad = n_BAdl_mat,sample_ids = which(sample_classes == 3),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[3],0.5),border = NA)
legend("topright",bty = "n",
       legend = c("Treatment 2","Treatment 1", "Control"),
       col = line_cols, lwd = 3, title = "Sub-adult")



# Tadpole
matplot(x = seq(1:ncol(n_BTad_mat)), y = t(n_BTad_mat),
        type = "l", xlab = "Species rank", ylab = "Abundance",
        main = "Log-log rank-abundance \n Tadpole", lwd = 1, lty=1,
        col=c(rep(line_cols[1], nrow(BTad_T2)),
              rep(line_cols[2], nrow(BTad_T1)),
              rep(line_cols[3], nrow(BTad_Ctr))),
        log = "xy")
# Metamorphic
matplot(x = seq(1:ncol(n_BMet_mat)), y = t(n_BMet_mat),
        type = "l", xlab = "Species rank", ylab = "Abundance",
        main = "Log-log rank-abundance \n Metamorphic", lwd = 1, lty=1,
        col=c(rep(line_cols[1], nrow(BMet_T2)),
              rep(line_cols[2], nrow(BMet_T1)),
              rep(line_cols[3], nrow(BMet_Ctr))),
        log = "xy")
# Sub-adult
matplot(x = seq(1:ncol(n_BAdl_mat)), y = t(n_BAdl_mat),
        type = "l", xlab = "Species rank", ylab = "Abundance",
        main = "Log-log rank-abundance \n Sub-adult", lwd = 1, lty=1,
        col=c(rep(line_cols[1], nrow(BAdl_T2)),
              rep(line_cols[2], nrow(BAdl_T1)),
              rep(line_cols[3], nrow(BAdl_Ctr))),
        log = "xy")


##### MSD analysis #####
# MSD analysis

# Tadpole
d <- dist(x = BTad_mat,method = "manhattan")
mds <- cmdscale(d = d,k = 5,eig = TRUE)
plot(mds$points,xlab = "First coordinate",ylab = "Second coordinate",pch = 19,
     cex =1,col = line_cols[sample_classes],
     main = "MDS analysis \n Tadpole")
sample_classes <- c(rep(1, nrow(BTad_T2)),rep(2, nrow(BTad_T1)),
                    rep(3, nrow(BTad_Ctr)))
a <- representative_point(input = mds$points,ids = which(sample_classes == 1),
                          col = scales::alpha(line_cols[1],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
a <- representative_point(input = mds$points,ids = which(sample_classes == 2),
                          col = scales::alpha(line_cols[2],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
a <- representative_point(input = mds$points,ids = which(sample_classes == 3),
                          col = scales::alpha(line_cols[3],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
legend("topleft",bty = "n",legend = c("Treatment 2","Treatment 1","Control"),
       col = line_cols,pch = 19)

# Metamorphic
d <- dist(x = BMet_mat,method = "manhattan")
mds <- cmdscale(d = d,k = 5,eig = TRUE)
sample_classes <- c(rep(1, nrow(BMet_T2)),rep(2, nrow(BMet_T2)),
                    rep(3, nrow(BMet_Ctr)))
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
legend("topleft",bty = "n",legend = c("Treatment 2","Treatment 1","Control"),
       col = line_cols,pch = 19)

# Sub-adults
d <- dist(x = BAdl_mat,method = "manhattan")
mds <- cmdscale(d = d,k = 5,eig = TRUE)
sample_classes <- c(rep(1, nrow(BAdl_T2)),rep(2, nrow(BAdl_T1)),
                    rep(3, nrow(BAdl_Ctr)))
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
legend("topleft",bty = "n",legend = c("Treatment 2","Treatment 1","Control"),
       col = line_cols,pch = 19)

#### Centrality analysis ####

ml_BTad <- ml_TaxGroup(ml_BTad, tax_bacter, "Phylum", "Genus")
ml_BMet <- ml_TaxGroup(ml_BMet, tax_bacter, "Phylum", "Genus")
ml_BAdl <- ml_TaxGroup(ml_BAdl, tax_bacter, "Phylum", "Genus")

BTad_degree <- ctr_df(ml_BTad, c("Treatment 2", "Treatment 1", "Control"))
BMet_degree <- ctr_df(ml_BMet, c("Treatment 2", "Treatment 1", "Control"))
BAdl_degree <- ctr_df(ml_BAdl, c("Treatment 2", "Treatment 1", "Control"))

BTad_phyl_degree <- phyl_ctr_df(BTad_degree, c("Treatment 2", "Treatment 1", "Control"),
                                n_layers = 3)
BMet_phyl_degree <- phyl_ctr_df(BMet_degree, c("Treatment 2", "Treatment 1", "Control"),
                                n_layers = 3)
BAdl_phyl_degree <- phyl_ctr_df(BAdl_degree, c("Treatment 2", "Treatment 1", "Control"),
                                n_layers = 3)

library(gridExtra)
grid.arrange(
  phyl_ctr_plot(BTad_phyl_degree, "Treatment 2", "Treatment 2"),
  phyl_ctr_plot(BTad_phyl_degree, "Treatment 1", "Treatment 1"),
  phyl_ctr_plot(BTad_phyl_degree, "Control", "Control"), ncol =1
)

grid.arrange(
  phyl_ctr_plot(BMet_phyl_degree, "Treatment 2", "Treatment 2"),
  phyl_ctr_plot(BMet_phyl_degree, "Treatment 1", "Treatment 1"),
  phyl_ctr_plot(BMet_phyl_degree, "Control", "Control"), ncol =1
)

grid.arrange(
  phyl_ctr_plot(BAdl_phyl_degree, "Treatment 2", "Treatment 2"),
  phyl_ctr_plot(BAdl_phyl_degree, "Treatment 1", "Treatment 1"),
  phyl_ctr_plot(BAdl_phyl_degree, "Control", "Control"), ncol =1
)




##### Multiplex centrality #####

###### Degree ######

ml_BTad_degree <- ctr_ml(ml_BTad, "degree")
ml_BMet_degree <- ctr_ml(ml_BMet, "degree")
ml_BAdl_degree <- ctr_ml(ml_BAdl, "degree")

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

###### Betweenness ######

ml_BTad_bet <- ctr_ml(ml_BTad, "betweenness")
ml_BMet_bet <- ctr_ml(ml_BMet, "betweenness")
ml_BAdl_bet <- ctr_ml(ml_BAdl, "betweenness")

# Tadpole
lay <- layoutMultiplex(ml_BTad_bet, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_BTad_bet, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_BTad, ml_BTad_bet, 2),
                 node.colors=node_color_mat(ml_BTad_bet, "centrality"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

# Metamorphic
lay <- layoutMultiplex(ml_BMet_bet, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_BMet_bet, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_BMet, ml_BMet_bet, 2),
                 node.colors=node_color_mat(ml_BMet_bet, "centrality"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

# Sub-adult
lay <- layoutMultiplex(ml_BAdl_bet, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_BAdl_bet, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_BAdl, ml_BAdl_bet, 2),
                 node.colors=node_color_mat(ml_BAdl_bet, "centrality"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

###### Closeness ######

ml_BTad_cl <- ctr_ml(ml_BTad, "closeness")
ml_BMet_cl <- ctr_ml(ml_BMet, "closeness")
ml_BAdl_cl <- ctr_ml(ml_BAdl, "closeness")

# Tadpole
lay <- layoutMultiplex(ml_BTad_cl, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_BTad_bet, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_BTad, ml_BTad_cl, 2),
                 node.colors=node_color_mat(ml_BTad_cl, "centrality"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

# Metamorphic
lay <- layoutMultiplex(ml_BMet_bet, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_BMet_bet, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_BMet, ml_BMet_bet, 2),
                 node.colors=node_color_mat(ml_BMet_bet, "centrality"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

# Sub-adult
lay <- layoutMultiplex(ml_BAdl_bet, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_BAdl_bet, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_BAdl, ml_BAdl_bet, 2),
                 node.colors=node_color_mat(ml_BAdl_bet, "centrality"),
                 edge.colors="#838B8B",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

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

#### Abundance vs degree ####

layer_names<-c("Control", "Treatment 1", "Treatment 2")
layer_colors<-c("green", "orange", "red")

cor_degree_abs(ml_BTad, abs_BTad, layer_names, layer_colors, "Tadpole")
cor_degree_abs(ml_BMet, abs_BMet, layer_names, layer_colors, "Metamorphic")
cor_degree_abs(ml_BAdl, abs_BAdl, layer_names, layer_colors, "Sub-adult")



#### Topological properties ####

treatments <- c(20, 25, 29)
frog_ml_properties<-rbind(ml_properties(ml_BTad, treatments),
                          ml_properties(ml_BMet, treatments),
                          ml_properties(ml_BAdl, treatments))

Stage = c(rep("Tadpole", length(ml_BTad)),
          rep("Metamorphic",length(ml_BTad)),
          rep("Sub-adult", length(ml_BTad)))

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

kable(slope_R2(ml_properties(ml_BTad, treatments))) %>% kable_styling()
kable(slope_R2(ml_properties(ml_BMet, treatments))) %>% kable_styling()
kable(slope_R2(ml_properties(ml_BAdl, treatments))) %>% kable_styling()

#### Dysbiosis analysis ####

# Load required packages
library(dysbiosisR)
library(ggplot2)
library(microbiome)
library(dplyr)
library(phyloseq)

# Transform the data into a phyloseq object
tax_bacter<-as.matrix(bacterias[,t1:t2])
otus_bacter<-as.matrix(bacterias[,-(t1:t2)])
otus_bacter<-otus_bacter[,-1]
# OTUs IDs as row names and sample IDs as column names
ID_otus <- bacterias$id
otus_bacter<-otus_bacter[,-1]
samples_bacter<-colnames(otus_bacter)
otus_bacter<-matrix(as.numeric(otus_bacter),
                    ncol = length(samples_bacter),
                    nrow = length(ID_otus))
colnames(otus_bacter)<-samples_bacter
rownames(otus_bacter)<-ID_otus
# OTUs IDs as row names in the taxa table
rownames(tax_bacter)<-ID_otus
# Phyloseq objects
otus_bacter<-otu_table(otus_bacter, taxa_are_rows = TRUE)
tax_bacter<-tax_table(tax_bacter)
physeq = phyloseq(otus_bacter, tax_bacter)
# Sample metadata
rownames(meta_bacterias)<-meta_bacterias[,1]
meta_bacterias<-meta_bacterias[,-1]
sample_data<-sample_data(meta_bacterias)
# Merge of phyloseq objects
physeq2<-merge_phyloseq(physeq, sample_data)

# Dysbiosis analysis

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

# Test for outliers’ detection that accounts for the wide range of
# microbiome phenotypes observed in a typical set of healthy individuals
# and for intra-individual temporal variation.
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


T2<-degree(ml_BAdl_sp[[3]])[-which(degree(ml_BAdl_sp[[3]])==0)]
T1<-degree(ml_BAdl_sp[[2]])[-which(degree(ml_BAdl_sp[[2]])==0)]
Ctr<-degree(ml_BAdl_sp[[1]])[-which(degree(ml_BAdl_sp[[1]])==0)]

T2<-degree(BTad_T2Net)[-which(degree(BTad_T2Net)==0)]
T1<-degree(BTad_T1Net)[-which(degree(BTad_T1Net)==0)]
Ctr<-degree(BTad_CtrNet)[-which(degree(BTad_CtrNet)==0)]

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
  title = "Sub-adults"
)

T2<-degree(BMet_T2Net)[-which(degree(BMet_T2Net)==0)]
T1<-degree(BMet_T1Net)[-which(degree(BMet_T1Net)==0)]
Ctr<-degree(BMet_CtrNet)[-which(degree(BMet_CtrNet)==0)]

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

T2<-degree(BAdl_T2Net)[-which(degree(BAdl_T2Net)==0)]
T1<-degree(BAdl_T1Net)[-which(degree(BAdl_T1Net)==0)]
Ctr<-degree(BAdl_CtrNet)[-which(degree(BAdl_CtrNet)==0)]

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


g<-erdos.renyi.game(100, 0.4)




