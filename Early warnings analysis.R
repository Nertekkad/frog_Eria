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

[which(colnames(bacterias) %in% adlID)]


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



mi_vector <- c(rnorm(2000, mean = 0, sd = 30), rnorm(2000, mean = 60, sd = 200))
diferencias <- diff(mi_vector)
plot(mi_vector[-length(mi_vector)], diferencias, type = "p", pch = 16,
     col = "blue", xlab = "Observaciones", ylab = "Diferencias")

















library(ggplot2)
library(dplyr)

# Simulación de datos (reemplaza esto con tus datos reales)
set.seed(123)
rr_intervals <- rnorm(500, mean = 1000, sd = 50)

# Calcula los intervalos RR consecutivos
rr_diff <- diff(rr_intervals)

# Crea un data frame con los intervalos RR y sus diferencias
poincare_data <- data.frame(RRn = rr_intervals[-length(rr_intervals)],
                            RRn1 = rr_intervals[-1],
                            RR_diff = rr_diff)

# Gráfico de Poincaré
ggplot(poincare_data, aes(x = RRn, y = RRn1)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(x = "RRn (ms)", y = "RRn+1 (ms)",
       title = "Diagrama de Poincaré",
       subtitle = "SD1 y SD2 se derivan de esta elipse") +
  theme_minimal()




