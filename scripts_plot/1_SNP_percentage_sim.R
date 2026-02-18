path <- getwd()
load(paste( path,"/simulation/Simulation_results/block_L_accuracy_128_sd1.RData", 
            sep=""))

colors= c(  "#FFC20A","#D41159" )


df_simu <- do.call(rbind, res)

colnames(df_simu) <- c("n_cs_nps",
                       "n_effect_nps",
                       "n_false_effect_nps",
                       "mean_purity_nps",
                       "cs_size_nps",
                       "n_cs_ps",
                       "n_effect_ps",
                       "n_false_effect_ps",
                       "mean_purity_ps",
                       "cs_size_ps",
                       "Number_effect",
                       "reg_sim")
df_simu <- as.data.frame(df_simu)
df_simu <- df_simu[which(df_simu$Number_effect %in% c(1,2,4,8,12,16)),]
df_simu <- df_simu[-which(df_simu$cs_size_ps>10),]
library(dplyr)
df_simu$power_nps <- df_simu$n_effect_nps/df_simu$Number_effect
df_simu$power_ps <- df_simu$n_effect_ps/df_simu$Number_effect
df_simu$t1_nps <- 1- df_simu$n_false_effect_nps/(df_simu$n_effect_nps+df_simu$n_false_effect_nps)
df_simu$t1_ps <- 1- df_simu$n_false_effect_ps/(df_simu$n_effect_ps+df_simu$n_false_effect_ps)



dft= df_simu
df1= df_simu[which(df_simu$cs_size_ps == 1),]
# 20% -----

### For PVE=20% ------
path <- getwd()
load(paste( path,"/simulation/Simulation_results/block_L_accuracy_128_sd2.RData", 
            sep=""))


df_simu <- do.call(rbind, res)

colnames(df_simu) <- c("n_cs_nps",
                       "n_effect_nps",
                       "n_false_effect_nps",
                       "mean_purity_nps",
                       "cs_size_nps",
                       "n_cs_ps",
                       "n_effect_ps",
                       "n_false_effect_ps",
                       "mean_purity_ps",
                       "cs_size_ps",
                       "Number_effect",
                       "reg_sim")
df_simu <- as.data.frame(df_simu)
df_simu <- df_simu[which(df_simu$Number_effect %in% c(1,2,4,8,12,16)),]
df_simu <- df_simu[-which(df_simu$cs_size_ps>10),]
library(dplyr)
df_simu$power_nps <- df_simu$n_effect_nps/df_simu$Number_effect
df_simu$power_ps <- df_simu$n_effect_ps/df_simu$Number_effect
df_simu$t1_nps <- 1- df_simu$n_false_effect_nps/(df_simu$n_effect_nps+df_simu$n_false_effect_nps)
df_simu$t1_ps <- 1- df_simu$n_false_effect_ps/(df_simu$n_effect_ps+df_simu$n_false_effect_ps)
df0=df_simu[which(df_simu$cs_size_ps == 1),]


df1= rbind(df1,df0)


### For PVE=30% ----- 
path <- getwd()
load(paste( path,"/simulation/Simulation_results/block_L_accuracy_128_sd3.RData", 
            sep=""))



df_simu <- do.call(rbind, res)

colnames(df_simu) <- c("n_cs_nps",
                       "n_effect_nps",
                       "n_false_effect_nps",
                       "mean_purity_nps",
                       "cs_size_nps",
                       "n_cs_ps",
                       "n_effect_ps",
                       "n_false_effect_ps",
                       "mean_purity_ps",
                       "cs_size_ps",
                       "Number_effect",
                       "reg_sim")
df_simu <- as.data.frame(df_simu)
df_simu <- df_simu[which(df_simu$Number_effect %in% c(1,2,4,8,12,16)),]
df_simu <- df_simu[-which(df_simu$cs_size_ps>10),]
library(dplyr)
df_simu$power_nps <- df_simu$n_effect_nps/df_simu$Number_effect
df_simu$power_ps <- df_simu$n_effect_ps/df_simu$Number_effect
df_simu$t1_nps <- 1- df_simu$n_false_effect_nps/(df_simu$n_effect_nps+df_simu$n_false_effect_nps)
df_simu$t1_ps <- 1- df_simu$n_false_effect_ps/(df_simu$n_effect_ps+df_simu$n_false_effect_ps)

df0=df_simu[which(df_simu$cs_size_ps == 1),]


df1= rbind(df1,df0)
dft= rbind(dft,df_simu)
### For PVE=30% ----- 

path <- getwd()
load(paste( path,"/simulation/Simulation_results/block_L_accuracy_128_sd4.RData", 
            sep=""))


df_simu <- do.call(rbind, res)
colnames(df_simu) <- c("n_cs_nps",
                       "n_effect_nps",
                       "n_false_effect_nps",
                       "mean_purity_nps",
                       "cs_size_nps",
                       "n_cs_ps",
                       "n_effect_ps",
                       "n_false_effect_ps",
                       "mean_purity_ps",
                       "cs_size_ps",
                       "Number_effect",
                       "reg_sim")
df_simu <- as.data.frame(df_simu)
df_simu <- df_simu[which(df_simu$Number_effect %in% c(1,2,4,8,12,16)),]
df_simu <- df_simu[-which(df_simu$cs_size_ps>10),]
library(dplyr)
df_simu$power_nps <- df_simu$n_effect_nps/df_simu$Number_effect
df_simu$power_ps <- df_simu$n_effect_ps/df_simu$Number_effect
df_simu$t1_nps <- 1- df_simu$n_false_effect_nps/(df_simu$n_effect_nps+df_simu$n_false_effect_nps)
df_simu$t1_ps <- 1- df_simu$n_false_effect_ps/(df_simu$n_effect_ps+df_simu$n_false_effect_ps)


df0=df_simu[which(df_simu$cs_size_ps == 1),]


dft= rbind(dft,df_simu)
df1= rbind(df1,df0)

## res ----

library(dplR)
### For PVE=10% ----

path <- getwd()
load(paste( path,"/simulation/Simulation_results/block_L_accuracy_128_sd1.RData", 
            sep=""))

colors= c(  "#FFC20A","#D41159" )


df_simu <- do.call(rbind, res)

colnames(df_simu) <- c("n_cs_nps",
                       "n_effect_nps",
                       "n_false_effect_nps",
                       "mean_purity_nps",
                       "cs_size_nps",
                       "n_cs_ps",
                       "n_effect_ps",
                       "n_false_effect_ps",
                       "mean_purity_ps",
                       "cs_size_ps",
                       "Number_effect",
                       "reg_sim")
df_simu <- as.data.frame(df_simu)
df_simu <- df_simu[which(df_simu$Number_effect %in% c(1,2,4,8,12,16)),]
df_simu <- df_simu[-which(df_simu$cs_size_ps>10),]
library(dplyr)
df_simu$power_nps <- df_simu$n_effect_nps/df_simu$Number_effect
df_simu$power_ps <- df_simu$n_effect_ps/df_simu$Number_effect
df_simu$t1_nps <- 1- df_simu$n_false_effect_nps/(df_simu$n_effect_nps+df_simu$n_false_effect_nps)
df_simu$t1_ps <- 1- df_simu$n_false_effect_ps/(df_simu$n_effect_ps+df_simu$n_false_effect_ps)

df0=df_simu[which(df_simu$cs_size_ps == 1),]


dft= rbind(dft,df_simu)
df1= rbind(df1,df0)

library(dplR)
### For PVE=20% ----

path <- getwd()
load(paste( path,"/simulation/Simulation_results/block_L_accuracy_128_sd2.RData", 
            sep=""))

colors= c(  "#FFC20A","#D41159" )


df_simu <- do.call(rbind, res)

colnames(df_simu) <- c("n_cs_nps",
                       "n_effect_nps",
                       "n_false_effect_nps",
                       "mean_purity_nps",
                       "cs_size_nps",
                       "n_cs_ps",
                       "n_effect_ps",
                       "n_false_effect_ps",
                       "mean_purity_ps",
                       "cs_size_ps",
                       "Number_effect",
                       "reg_sim")
df_simu <- as.data.frame(df_simu)
df_simu <- df_simu[which(df_simu$Number_effect %in% c(1,2,4,8,12,16)),]
df_simu <- df_simu[-which(df_simu$cs_size_ps>10),]
library(dplyr)
df_simu$power_nps <- df_simu$n_effect_nps/df_simu$Number_effect
df_simu$power_ps <- df_simu$n_effect_ps/df_simu$Number_effect
df_simu$t1_nps <- 1- df_simu$n_false_effect_nps/(df_simu$n_effect_nps+df_simu$n_false_effect_nps)
df_simu$t1_ps <- 1- df_simu$n_false_effect_ps/(df_simu$n_effect_ps+df_simu$n_false_effect_ps)

df0=df_simu[which(df_simu$cs_size_ps == 1),]

dft= rbind(dft,df_simu)

df1= rbind(df1,df0)

library(dplR)
### For PVE=30% ----

path <- getwd()
load(paste( path,"/simulation/Simulation_results/block_L_accuracy_128_sd3.RData", 
            sep=""))

colors= c(  "#FFC20A","#D41159" )


df_simu <- do.call(rbind, res)

colnames(df_simu) <- c("n_cs_nps",
                       "n_effect_nps",
                       "n_false_effect_nps",
                       "mean_purity_nps",
                       "cs_size_nps",
                       "n_cs_ps",
                       "n_effect_ps",
                       "n_false_effect_ps",
                       "mean_purity_ps",
                       "cs_size_ps",
                       "Number_effect",
                       "reg_sim")
df_simu <- as.data.frame(df_simu)
df_simu <- df_simu[which(df_simu$Number_effect %in% c(1,2,4,8,12,16)),]
df_simu <- df_simu[-which(df_simu$cs_size_ps>10),]
library(dplyr)
df_simu$power_nps <- df_simu$n_effect_nps/df_simu$Number_effect
df_simu$power_ps <- df_simu$n_effect_ps/df_simu$Number_effect
df_simu$t1_nps <- 1- df_simu$n_false_effect_nps/(df_simu$n_effect_nps+df_simu$n_false_effect_nps)
df_simu$t1_ps <- 1- df_simu$n_false_effect_ps/(df_simu$n_effect_ps+df_simu$n_false_effect_ps)

df0=df_simu[which(df_simu$cs_size_ps == 1),]


dft= rbind(dft,df_simu)
df1= rbind(df1,df0)

library(dplR)
### For PVE=40% ----

path <- getwd()
load(paste( path,"/simulation/Simulation_results/block_L_accuracy_128_sd4.RData", 
            sep=""))

colors= c(  "#FFC20A","#D41159" )


df_simu <- do.call(rbind, res)

colnames(df_simu) <- c("n_cs_nps",
                       "n_effect_nps",
                       "n_false_effect_nps",
                       "mean_purity_nps",
                       "cs_size_nps",
                       "n_cs_ps",
                       "n_effect_ps",
                       "n_false_effect_ps",
                       "mean_purity_ps",
                       "cs_size_ps",
                       "Number_effect",
                       "reg_sim")
df_simu <- as.data.frame(df_simu)
df_simu <- df_simu[which(df_simu$Number_effect %in% c(1,2,4,8,12,16)),]
df_simu <- df_simu[-which(df_simu$cs_size_ps>10),]
library(dplyr)
df_simu$power_nps <- df_simu$n_effect_nps/df_simu$Number_effect
df_simu$power_ps <- df_simu$n_effect_ps/df_simu$Number_effect
df_simu$t1_nps <- 1- df_simu$n_false_effect_nps/(df_simu$n_effect_nps+df_simu$n_false_effect_nps)
df_simu$t1_ps <- 1- df_simu$n_false_effect_ps/(df_simu$n_effect_ps+df_simu$n_false_effect_ps)

df0=df_simu[which(df_simu$cs_size_ps == 1),]


dft= rbind(dft,df_simu)
df1= rbind(df1,df0)

df1

dft
