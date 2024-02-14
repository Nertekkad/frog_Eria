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
FTad_T1Net<-net_inference(FTad_T1, "sparcc")
FTad_T2Net<-net_inference(FTad_T2, "sparcc")
FTad_CtrNet<-net_inference(FTad_Ctr, "sparcc")
# Metamorphic
FMet_T1Net<-net_inference(FMet_T1, "sparcc")
FMet_T2Net<-net_inference(FMet_T2, "sparcc")
FMet_CtrNet<-net_inference(FMet_Ctr, "sparcc")
# Sub-adult
FAdl_T1Net<-net_inference(FAdl_T1, "sparcc")
FAdl_T2Net<-net_inference(FAdl_T2, "sparcc")
FAdl_CtrNet<-net_inference(FAdl_Ctr, "sparcc")
