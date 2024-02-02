#### Data of frog's microbiota under different temperatures ####

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

unq<-unique(tax_bacter[,"Phylum"])
unq<-unq[-c(which(unq == "uncultured"), which(is.na(unq)))]
colors <- sample(rainbow(100), length(unq))

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

plot(BTad_T1Net, vertex.label.color="black",
     vertex.color = vertex.attributes(BTad_T1Net)$color, vertex.label.cex=.5,
     vertex.label.dist=1,layout=layout_with_kk, vertex.size = 5,
     main = "Tadpole under treatment 1")
legend(x=-2.4, y=0.7, unq, title = "Insect", pch=21, pt.bg=colors, pt.cex=1.3, cex=.8, bty="n", ncol=1)


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
n_BTad_mat<-n_BTad_mat[,-which(colSums(n_BTad_mat)==0)]

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

# Tadpole
d <- dist(x = BTad_mat,method = "manhattan")
mds <- cmdscale(d = d,k = 5,eig = TRUE)
plot(mds$points,xlab = "First coordinate",ylab = "Second coordinate",pch = 19,
     cex =1,col = line_cols[sample_classes],
     main = "MDS plot with representative points \n of each group and error bars")
sample_classes <- c(rep(1, length(T2_s)),rep(2, length(T1_s)),
                    rep(3, length(TC_s)))
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

# Metamorphic
d <- dist(x = n_BMet_mat,method = "manhattan")
mds <- cmdscale(d = d,k = 5,eig = TRUE)
sample_classes <- c(rep(1, length(M2_s)),rep(2, length(M1_s)),
                    rep(3, length(MC_s)))
plot(mds$points,xlab = "First coordinate",ylab = "Second coordinate",pch = 19,
     cex =1,col = line_cols[sample_classes],
     main = "MDS analysis \n ")
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

# Sub-adults
d <- dist(x = n_BAdl_mat,method = "manhattan")
mds <- cmdscale(d = d,k = 5,eig = TRUE)
sample_classes <- c(rep(1, length(A2_s)),rep(2, length(A1_s)),
                    rep(3, length(AC_s)))
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
ml_BTad <- ctr_ml(ml_BTad, "degree")
ml_BMet <- ctr_ml(ml_BMet, "degree")
ml_BAdl <- ctr_ml(ml_BAdl, "degree")




degree_df <- data.frame(Species = vertex.attributes(BTad_CtrNet)$name,
                        color = vertex.attributes(BTad_CtrNet)$color,
                        Phylum = vertex.attributes(BTad_CtrNet)$Taxon,
                        Ctr_degree = degree(BTad_CtrNet),
                        T1_degree = degree(BTad_T1Net),
                        T2_degree = degree(BTad_T2Net))

degree_df<-degree_df[-which(rowSums(degree_df[,c(4,5,6)])/
                              mean(rowSums(degree_df[,c(4,5,6)]))<0.9),]

library(gridExtra)
p1<-ggplot(data=degree_df, aes(x=Species, y=Ctr_degree, fill=Phylum)) +
  geom_bar(stat="identity")
p2<-ggplot(data=degree_df, aes(x=Species, y=T1_degree, fill=Phylum)) +
  geom_bar(stat="identity")
p3<-ggplot(data=degree_df, aes(x=Species, y=T2_degree, fill=Phylum)) +
  geom_bar(stat="identity")

grid.arrange(p1 + coord_flip() +
               theme(axis.text.y = element_text(size = 0.01),
                     legend.position = "none") + ylab("Control"),
             p2 + coord_flip() +
               theme(axis.text.y = element_text(size = 0.01),
                     legend.position = "none") + ylab("Treatment 1"),
             p3 + coord_flip() +
               theme(axis.text.y = element_text(size = 0.01),
                     legend.text = element_text(size = 6),
                     legend.position="bottom", legend.box = "horizontal") +
               ylab("Treatment 2"),
             ncol=3)