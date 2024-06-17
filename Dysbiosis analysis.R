#### Dysbiosis analysis ####

#### Data of frog's bacteriome under different temperatures ####

setwd("~/Frog's data")
# Data loading
library(readxl)
meta_bacterias<-read.csv("Metadata_Bacterias.csv")
bacterias<-read_excel("Bacterias.xlsx")

# Load required packages
library(dysbiosisR)
library(ggplot2)
library(microbiome)
library(dplyr)
library(phyloseq)


#### Bacteria ####

##### Tadpoles ######
Tadmet<-meta_bacterias[which(meta_bacterias$Life_stage %in% "Tadpole"),]
rownames(Tadmet)<-Tadmet[,1]
Tadotu<-bacterias[,which(colnames(bacterias) %in% rownames(Tadmet))]
rownames(Tadotu)<-bacterias$id
Tadmet_phy<-sample_data(Tadmet)
Tadotu_phy<-otu_table(Tadotu, T)
physeq2<-merge_phyloseq(Tadotu_phy, Tadmet_phy)



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


library(MicrobiotaProcess)
library(UpSetR)

upsetda2 <- get_upset(obj=physeq2, factorNames="Treatment")
upset(upsetda2, sets=unique(as.vector(sample_data(physeq2)$Treatment)), 
      sets.bar.color = "#56B4E9",
      order.by = "freq", 
      empty.intersections = "on")


##### Metamorphic #####

Metmet<-meta_bacterias[which(meta_bacterias$Life_stage %in% "Metamorphic"),]
rownames(Metmet)<-Metmet[,1]
Metotu<-bacterias[,which(colnames(bacterias) %in% rownames(Metmet))]
rownames(Metotu)<-bacterias$id
Metmet_phy<-sample_data(Metmet)
Metotu_phy<-otu_table(Metotu, T)
physeq2<-merge_phyloseq(Metotu_phy, Metmet_phy)




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

library(MicrobiotaProcess)
library(UpSetR)

upsetda2 <- get_upset(obj=physeq2, factorNames="Treatment")
upset(upsetda2, sets=unique(as.vector(sample_data(physeq2)$Treatment)), 
      sets.bar.color = "#56B4E9",
      order.by = "freq", 
      empty.intersections = "on")


##### Sub-adults #####

Adlmet<-meta_bacterias[which(meta_bacterias$Life_stage %in% "Sub-adult"),]
rownames(Adlmet)<-Adlmet[,1]
Adlotu<-bacterias[,which(colnames(bacterias) %in% rownames(Adlmet))]
rownames(Adlotu)<-bacterias$id
Adlmet_phy<-sample_data(Adlmet)
Adlotu_phy<-otu_table(Adlotu, T)
physeq2<-merge_phyloseq(Adlotu_phy, Adlmet_phy)


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

library(MicrobiotaProcess)
library(UpSetR)

upsetda2 <- get_upset(obj=physeq2, factorNames="Treatment")
upset(upsetda2, sets=unique(as.vector(sample_data(physeq2)$Treatment)), 
      sets.bar.color = "#56B4E9",
      order.by = "freq", 
      empty.intersections = "on")

#### Fungus ####

# Data loading
meta_hongos<-read.csv("Metadata_Hongos.csv")
hongos<-read_excel("Hongos.xlsx")

##### Tadpoles ######
Tadmet<-meta_hongos[which(meta_hongos$Life_stage %in% "Tadpole"),]
rownames(Tadmet)<-Tadmet[,1]
Tadotu<-hongos[,which(colnames(hongos) %in% rownames(Tadmet))]
rownames(Tadotu)<-hongos$id
Tadmet_phy<-sample_data(Tadmet)
Tadotu_phy<-otu_table(Tadotu, T)
physeq2<-merge_phyloseq(Tadotu_phy, Tadmet_phy)

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


##### Metamorphic ######
Tadmet<-meta_hongos[which(meta_hongos$Life_stage %in% "Metamorphic"),]
rownames(Tadmet)<-Tadmet[,1]
Tadotu<-hongos[,which(colnames(hongos) %in% rownames(Tadmet))]
rownames(Tadotu)<-hongos$id
Tadmet_phy<-sample_data(Tadmet)
Tadotu_phy<-otu_table(Tadotu, T)
physeq2<-merge_phyloseq(Tadotu_phy, Tadmet_phy)

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


##### Sub-adults ######
Tadmet<-meta_hongos[which(meta_hongos$Life_stage %in% "Sub-adult"),]
rownames(Tadmet)<-Tadmet[,1]
Tadotu<-hongos[,which(colnames(hongos) %in% rownames(Tadmet))]
rownames(Tadotu)<-hongos$id
Tadmet_phy<-sample_data(Tadmet)
Tadotu_phy<-otu_table(Tadotu, T)
physeq2<-merge_phyloseq(Tadotu_phy, Tadmet_phy)

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

