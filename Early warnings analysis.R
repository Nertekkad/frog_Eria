library(earlywarnings)
library(EWS)
library(EWSmethods)
library(codyn)
library(vegan)


setwd("~/Frog's data")
# Data loading
library(readxl)
meta_bacterias<-read.csv("Metadata_Bacterias.csv")
bacterias<-read_excel("Bacterias.xlsx")

tadID<-meta_bacterias$sample.id[which(meta_bacterias$Life_stage %in% "Tadpole")]
metID<-meta_bacterias$sample.id[which(meta_bacterias$Life_stage %in% "Metamorphic")]
adlID<-meta_bacterias$sample.id[which(meta_bacterias$Life_stage %in% "Sub-adult")]

#### Richness plots ####

richness<-function(ab_table){
  richness<-c()
  for(i in 1:nrow(ab_table)){
    richness[i]<-ab_table[i,][which(ab_table[i,]!=0)]
  }
  return(richness)
}

richness(BTad_Ctr)
richness(BTad_T1)
richness(BTad_T2)




# Tadpole
df <- data.frame(x = 1:length(richness(BTad_Ctr)), y = richness(BTad_Ctr))
ggplot(df, aes(x, y)) +
  geom_line(size = 1.5, color = "darkblue") +
  labs(x = "Samples", y = "Richness", title = "Control")

df <- data.frame(x = 1:length(richness(BTad_T1)), y = richness(BTad_T1))
ggplot(df, aes(x, y)) +
  geom_line(size = 1.5, color = "darkblue") +
  labs(x = "Samples", y = "Richness", title = "Treatment 1")

df <- data.frame(x = 1:length(richness(BTad_T2)), y = richness(BTad_T2))
ggplot(df, aes(x, y)) +
  geom_line(size = 1.5, color = "darkblue") +
  labs(x = "Samples", y = "Richness", title = "Treatment 2")

# Metamorphic
df <- data.frame(x = 1:length(richness(BMet_Ctr)), y = richness(BMet_Ctr))
ggplot(df, aes(x, y)) +
  geom_line(size = 1.5, color = "darkblue") +
  labs(x = "Samples", y = "Richness", title = "Control")

df <- data.frame(x = 1:length(richness(BMet_T1)), y = richness(BMet_T1))
ggplot(df, aes(x, y)) +
  geom_line(size = 1.5, color = "darkblue") +
  labs(x = "Samples", y = "Richness", title = "Treatment 1")

df <- data.frame(x = 1:length(richness(BMet_T2)), y = richness(BMet_T2))
ggplot(df, aes(x, y)) +
  geom_line(size = 1.5, color = "darkblue") +
  labs(x = "Samples", y = "Richness", title = "Treatment 2")

# Sub-adult
df <- data.frame(x = 1:length(richness(BAdl_Ctr)), y = richness(BAdl_Ctr))
ggplot(df, aes(x, y)) +
  geom_line(size = 1.5, color = "darkblue") +
  labs(x = "Samples", y = "Richness", title = "Control")

df <- data.frame(x = 1:length(richness(BAdl_T1)), y = richness(BAdl_T1))
ggplot(df, aes(x, y)) +
  geom_line(size = 1.5, color = "darkblue") +
  labs(x = "Samples", y = "Richness", title = "Treatment 1")

df <- data.frame(x = 1:length(richness(BAdl_T2)), y = richness(BAdl_T2))
ggplot(df, aes(x, y)) +
  geom_line(size = 1.5, color = "darkblue") +
  labs(x = "Samples", y = "Richness", title = "Treatment 2")

# Tadpole
df_violinplot<-data.frame(
  Treatment = c(rep("Treatment 2", length(richness(BTad_T2))),
                rep("Treatment 1", length(richness(BTad_T1))),
                rep("Control", length(richness(BTad_Ctr)))),
  Data = c(richness(BTad_T2), richness(BTad_T1), richness(BTad_Ctr))
)

library(ggstatsplot)
ggbetweenstats(
  data  = df_violinplot,
  x     = Treatment,
  y     = Data,
  title = "Tadpoles"
)

# Metamorphic
df_violinplot<-data.frame(
  Treatment = c(rep("Treatment 2", length(richness(BMet_T2))),
                rep("Treatment 1", length(richness(BMet_T1))),
                rep("Control", length(richness(BMet_Ctr)))),
  Data = c(richness(BMet_T2), richness(BMet_T1), richness(BMet_Ctr))
)

library(ggstatsplot)
ggbetweenstats(
  data  = df_violinplot,
  x     = Treatment,
  y     = Data,
  title = "Metamorphic"
)

# Sub-adult
df_violinplot<-data.frame(
  Treatment = c(rep("Treatment 2", length(richness(BAdl_T2))),
                rep("Treatment 1", length(richness(BAdl_T1))),
                rep("Control", length(richness(BAdl_Ctr)))),
  Data = c(richness(BAdl_T2), richness(BAdl_T1), richness(BAdl_Ctr))
)

library(ggstatsplot)
ggbetweenstats(
  data  = df_violinplot,
  x     = Treatment,
  y     = Data,
  title = "Sub-adult"
)

#### Early warnings metrics ####

# Shannon diversity
Div_Shannon <- function(abundancias_ab){
  abs_rel <- abundancias_ab/sum(abundancias_ab)
  Shannon <- -sum(abs_rel*log(abs_rel))
  return(Shannon)
}

# Simpson dominance
Dom_Simpson <- function(abundancias_ab){
  abs_rel <- abundancias_ab/sum(abundancias_ab)
  Simpson <- sum(abs_rel^2)
  return(Simpson)
}

# Pielou evenness
Eq_Pielou <- function(abundancias_ab){
  abs_rel <- abundancias_ab/sum(abundancias_ab)
  Shannon <- -sum(abs_rel*log(abs_rel))
  Pielou <- Shannon/log(length(abundancias_ab))
  return(Pielou)
}

# Diversity function for abundance tables
ab_table_div<-function(ab_table, diversity_type){
  require(vegan)
  if(diversity_type == "shannon"){
    div_table<-c()
    for(i in 1:ncol(ab_table)){
      div_table[i]<-diversity(ab_table[, i])
    }
  } else
    if(diversity_type == "simpson"){
      div_table<-c()
      for(i in 1:ncol(ab_table)){
        div_table[i]<-Dom_Simpson(ab_table[, i])
      }
    } else
      if(diversity_type == "pielou"){
        div_table<-c()
        for(i in 1:ncol(ab_table)){
          S <- length(ab_table[, i])
          div_table[i] <- diversity(ab_table[, i])/log(S)
        }
      } else
        if(diversity_type == "ginisimpson"){
          div_table<-c()
          for(i in 1:ncol(ab_table)){
            div_table[i]<-1-Dom_Simpson(ab_table[, i])
          }
        }
  return(div_table)
}


BTad<-rbind(BTad_Ctr, BTad_T1, BTad_T2)
BMet<-rbind(BMet_Ctr, BMet_T1, BMet_T2)
BAdl<-rbind(BAdl_Ctr, BAdl_T1, BAdl_T2)

# Tadpole

Tad_data <- data.frame(time = seq(1:length(ab_table_div(t(BTad), "shannon"))),
                       abundance = ab_table_div(t(BTad), "shannon")) #dummy skylark dataset

ews_metrics <- c("SD","ar1","skew") #the early warning signal metrics we wish to compute

roll_ews <- uniEWS(data = Tad_data, metrics =  ews_metrics, method = "rolling", winsize = 50) #lets use a rolling window approach

roll_ews$EWS$cor
plot(roll_ews,  y_lab = "Abundances")

exp_ews <- uniEWS(data = Tad_data, metrics =  ews_metrics, method = "expanding",
                  burn_in = 10, threshold = 2,  tail.direction = "one.tailed")
plot(exp_ews, y_lab = "Abundances")

# Metamorphic

Met_data <- data.frame(time = seq(1:length(ab_table_div(t(BMet), "shannon"))),
                       abundance = ab_table_div(t(BMet), "shannon")) #dummy skylark dataset

ews_metrics <- c("SD","ar1","skew") #the early warning signal metrics we wish to compute

roll_ews <- uniEWS(data = Met_data, metrics =  ews_metrics, method = "rolling", winsize = 50) #lets use a rolling window approach

roll_ews$EWS$cor

plot(roll_ews,  y_lab = "Abundances")

exp_ews <- uniEWS(data = Met_data, metrics =  ews_metrics, method = "expanding",
                  burn_in = 10, threshold = 2,  tail.direction = "one.tailed")
plot(exp_ews, y_lab = "Abundances")


# Sub-adult

Adl_data <- data.frame(time = seq(1:length(ab_table_div(t(BAdl), "shannon"))),
                       abundance = ab_table_div(t(BAdl), "shannon")) #dummy skylark dataset

ews_metrics <- c("SD","ar1","skew") #the early warning signal metrics we wish to compute

roll_ews <- uniEWS(data = Adl_data, metrics =  ews_metrics, method = "rolling", winsize = 50) #lets use a rolling window approach

roll_ews$EWS$cor

plot(roll_ews,  y_lab = "Abundances")

exp_ews <- uniEWS(data = Adl_data, metrics =  ews_metrics, method = "expanding",
                  burn_in = 10, threshold = 2,  tail.direction = "one.tailed")
plot(exp_ews, y_lab = "Abundances")




#### Auto-correlation ####

###### Separation of data by life stages type ######

# Identify frog´s life stages
life_stage<-unique(meta_bacterias$Life_stage)

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


poincare_plot<-function(ab_table, title){
  require(ggplot2)
  require(data.table)
  # Create vector
  data <- ab_table_div(ab_table, "shannon")
  data<-datos[-which(data == 0)]
  # Quantify the dots for Poincaré plot
  poincare_data <- data.frame(
    x = data[-length(data)], 
    y = data[-1]
  )
  # Create plot
  p <- ggplot(poincare_data, aes(x = x, y = y)) +
    geom_point() +
    theme_minimal()
  # Add ellipse
  p <- p + stat_ellipse(type = "t", level = 0.95) + labs(x = "RR_n", y = "RR_n+1", title = title)
  return(p)
}

pc_Tad_Ctr<-poincare_plot(BTad_Ctr, "Tadpole's control")
pc_Tad_T1<-poincare_plot(BTad_T1, "Tadpole's treatment 1")
pc_Tad_T2<-poincare_plot(BTad_T2, "Tadpole's treatment 2")

pc_Met_Ctr<-poincare_plot(BMet_Ctr, "Metamorphic's control")
pc_Met_T1<-poincare_plot(BMet_T1, "Metamorphoic's treatment 1")
pc_Met_T2<-poincare_plot(BMet_T2, "Metamorphic's treatment 2")

pc_Adl_Ctr<-poincare_plot(BAdl_Ctr, "Sub-adult's control")
pc_Adl_T1<-poincare_plot(BAdl_T1, "Sub-adult's treatment 1")
pc_Adl_T2<-poincare_plot(BAdl_T2, "Sub-adult's treatment 2")

library(gridExtra)
grid.arrange(pc_Tad_Ctr, pc_Tad_T1, pc_Tad_T2,
             pc_Met_Ctr, pc_Met_T1, pc_Met_T2,
             pc_Adl_Ctr, pc_Adl_T1, pc_Adl_T2,
             ncol = 3)
