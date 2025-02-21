rm(list=ls())

colors <- c("#D41159","#1A85FF","#40B0A6" )


### For PVE=10% ----
#### Calibration and  power  plot ----

#### Gaussian functional ---- 
path <- getwd()
load(paste( path,"/simulation/Simulation_results/check_L_accuracy_128_sd1.RData", 
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
#df_simu <- df_simu[which(df_simu$Number_effect %in% c(1,2,4,8,12,16)),]
df_simu <- df_simu[-which(df_simu$cs_size_ps>10),]
library(dplyr)
df_simu$power_nps <- df_simu$n_effect_nps/df_simu$Number_effect
df_simu$power_ps <- df_simu$n_effect_ps/df_simu$Number_effect
df_simu$t1_nps <- 1- df_simu$n_false_effect_nps/(df_simu$n_effect_nps+df_simu$n_false_effect_nps)
df_simu$t1_ps <- 1- df_simu$n_false_effect_ps/(df_simu$n_effect_ps+df_simu$n_false_effect_ps)



mean_power_nps <- rep( NA, max(df_simu$Number_effect) )
mean_power_ps <- rep( NA, max(df_simu$Number_effect) )
mean_T1_nps <- rep( NA,  max(df_simu$Number_effect) )
mean_T1_ps <- rep( NA,  max(df_simu$Number_effect) )

mean_cs_size_nps <- rep( NA,  max(df_simu$Number_effect) )
mean_cs_size_ps <- rep( NA,  max(df_simu$Number_effect) )

mean_purity_nps <- rep( NA,  max(df_simu$Number_effect) )
mean_purity_ps <- rep( NA,  max(df_simu$Number_effect) )
which_L <- rep( NA,  max(df_simu$Number_effect) )

h <- 1
for ( i in unique(df_simu$Number_effect))
{
  
  
  mean_power_nps[h] <- mean(df_simu$power_nps[which(df_simu$Number_effect ==i )] )
  mean_power_ps [h]<- mean(df_simu$power_ps[which(df_simu$Number_effect ==i)  ] )
  
  mean_T1_nps   [h] <- mean(df_simu$t1_nps[which(df_simu$Number_effect ==i  )] )
  
  mean_T1_ps    [h]<- mean(df_simu$t1_ps[which(df_simu$Number_effect ==i )] )
  
  
  mean_cs_size_nps[h] <- tbrm(df_simu$cs_size_nps[which(df_simu$Number_effect ==i )] )
  mean_cs_size_ps[h]  <- tbrm(df_simu$cs_size_ps[which(df_simu$Number_effect ==i )] )
  
  mean_purity_nps[h] <-  mean(df_simu$mean_purity_nps[which(df_simu$Number_effect ==i )] )
  mean_purity_ps[h]  <- mean(df_simu$mean_purity_ps[which(df_simu$Number_effect ==i )] )
  
  
  which_L       [h] <- i
  
  
  h <- h+1
  
  
  
  
}

final_df1 <- rbind(mean_power_nps ,    mean_power_ps,
                   mean_T1_nps ,mean_T1_ps ,
                   mean_cs_size_nps,
                   mean_cs_size_ps,
                   mean_purity_nps,
                   mean_purity_ps,
                   which_L  )

t_names <-rownames(final_df1)
final_df1 <- data.frame(t(final_df1))
colnames(final_df1) <- t_names



df_plot <- data.frame(power=c(final_df1$mean_power_ps, final_df1$mean_power_nps),
                      T1_error =c(final_df1$mean_T1_ps, final_df1$mean_T1_nps),
                      cs_size= c(final_df1$mean_cs_size_ps, final_df1$mean_cs_size_nps),
                      mean_purity =  c(final_df1$mean_purity_ps, final_df1$mean_purity_nps),
                      prior=as.factor(c(rep( "SPSP", nrow(final_df1)),rep( "ISP", nrow(final_df1)))),
                      L=  factor(rep( final_df1$which_L,2)))
df_plot <- df_plot[complete.cases(df_plot),]
sd_error_bin_up <- function(p,n_rep=100){
  c(1.96*sqrt(p*(1-p)/n_rep)+p)
}
sd_error_bin_low <- function(p,n_rep=100){
  c(-1.96*sqrt(p*(1-p)/n_rep)+p)
}

n_rep <- rep(c(table(df_simu$Number_effect)),each=2)

df_plot$pw_er_up <- sd_error_bin_up(df_plot[,1],n_rep)
df_plot$pw_er_low <- sd_error_bin_low(df_plot[,1],n_rep)

df_plot$t1_er_up <- sd_error_bin_up(df_plot[,2],n_rep)
df_plot$t1_er_low <- sd_error_bin_low(df_plot[,2],n_rep)

path <- getwd()
load(paste0(path,
            "/simulation/Simulation_script/script/additional_simualtion_for_fig3_panel_D_susie_pc_calibration_power/check_L_susie_128_sd1.RData"))





df_susie= data.frame(do.call( rbind, res))
colnames(df_susie)=c( "Number_effect", "n_cs", "n_effect", "n_false_effect")

df_simu=df_susie
df_simu$power  <- df_simu$n_effect /df_simu$Number_effect 
df_simu$t1  <- 1- df_simu$n_false_effect /(df_simu$n_effect +df_simu$n_false_effect )




mean_power_susie <- rep( NA, max(df_simu$Number_effect) ) 
mean_T1_susie <- rep( NA,  max(df_simu$Number_effect) )

which_L <- rep( NA,  max(df_simu$Number_effect) )

h <- 1
for ( i in unique(df_simu$Number_effect))
{
  
  
  mean_power_susie[h] <- mean(df_simu$power [which(df_simu$Number_effect ==i )] )
  
  mean_T1_susie   [h] <- mean(df_simu$t1 [which(df_simu$Number_effect ==i  )] )
  
  
  which_L       [h] <- i
  
  
  h <- h+1
  
  
  
  
}

mean_T1_susie
mean_power_susie


df= data.frame(power= mean_power_susie,
               T1_error= mean_T1_susie,
               
               L= unique(df_simu$Number_effect))


plot( df$L, df$Power)
plot( df$L, df$T1)


sd_error_bin_up <- function(p,n_rep=100){
  c(1.96*sqrt(p*(1-p)/n_rep)+p)
}
sd_error_bin_low <- function(p,n_rep=100){
  c(-1.96*sqrt(p*(1-p)/n_rep)+p)
}

n_rep <- rep(c(table(df_simu$Number_effect)),each=1)

df$pw_er_up <- sd_error_bin_up(df[,1],n_rep)
df$pw_er_low <- sd_error_bin_low(df[,1],n_rep)

df$t1_er_up <- sd_error_bin_up(df[,2],n_rep)
df$t1_er_low <- sd_error_bin_low(df[,2],n_rep)

df$L=as.factor(df$L)
df$prior=  rep( "SuSiE", nrow(df))


dfp =  rbind ( df_plot[,c("power",
                          "T1_error", 
                          "prior",
                          "L",
                          "pw_er_up",
                          "pw_er_low",
                          "t1_er_up",
                          "t1_er_low") ],
               df     [,c("power",
                          "T1_error", 
                          "prior",
                          "L",
                          "pw_er_up",
                          "pw_er_low",
                          "t1_er_up",
                          "t1_er_low") ])



library(ggplot2)

P11  <- ggplot(dfp, aes( x=L, y= power, col=prior))+
  geom_point(position=position_dodge(.9),size=2)+
  geom_errorbar(aes(ymin=pw_er_low, ymax=pw_er_up), width=.2,
                position=position_dodge(.9) )+
  ylim(c(0 ,1))+
  ylab("Power")+
  theme_linedraw()+
  scale_color_manual(values = colors) +
  scale_x_discrete(labels =c(1,"",3,"",5,"",7,"",9,"",11,"",13,"",15,"",17,"",19,""))

P11_t1 <- 
ggplot(dfp, aes( x=L, y= T1_error, col=prior))+
  
  geom_hline(yintercept = 0.95, color="black", linewidth=2)+
  geom_point(position=position_dodge(.9),size=2)+
  geom_errorbar(aes(ymin=t1_er_low, ymax=t1_er_up), width=.2,
                position=position_dodge(.9) )+
  ylim(c(0.7,1.04))+
  ylab("Coverage")+
  #ggtitle("Coverage PVE=10%") +
  theme_linedraw()+
  scale_color_manual(values = colors)+
  scale_x_discrete(labels =c(1,"",3,"",5,"",7,"",9,"",11,"",13,"",15,"",17,"",19,""))
P11_t1

## block ---- 

path <- getwd()
load(paste( path,"/simulation/Simulation_results/block_L_accuracy_128_sd1.RData", 
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
#df_simu <- df_simu[which(df_simu$Number_effect %in% c(1,2,4,8,12,16)),]
df_simu <- df_simu[-which(df_simu$cs_size_ps>10),]
library(dplyr)
df_simu$power_nps <- df_simu$n_effect_nps/df_simu$Number_effect
df_simu$power_ps <- df_simu$n_effect_ps/df_simu$Number_effect
df_simu$t1_nps <- 1- df_simu$n_false_effect_nps/(df_simu$n_effect_nps+df_simu$n_false_effect_nps)
df_simu$t1_ps <- 1- df_simu$n_false_effect_ps/(df_simu$n_effect_ps+df_simu$n_false_effect_ps)



mean_power_nps <- rep( NA, max(df_simu$Number_effect) )
mean_power_ps <- rep( NA, max(df_simu$Number_effect) )
mean_T1_nps <- rep( NA,  max(df_simu$Number_effect) )
mean_T1_ps <- rep( NA,  max(df_simu$Number_effect) )

mean_cs_size_nps <- rep( NA,  max(df_simu$Number_effect) )
mean_cs_size_ps <- rep( NA,  max(df_simu$Number_effect) )

mean_purity_nps <- rep( NA,  max(df_simu$Number_effect) )
mean_purity_ps <- rep( NA,  max(df_simu$Number_effect) )
which_L <- rep( NA,  max(df_simu$Number_effect) )

h <- 1
for ( i in unique(df_simu$Number_effect))
{
  
  
  mean_power_nps[h] <- mean(df_simu$power_nps[which(df_simu$Number_effect ==i )] )
  mean_power_ps [h]<- mean(df_simu$power_ps[which(df_simu$Number_effect ==i)  ] )
  
  mean_T1_nps   [h] <- mean(df_simu$t1_nps[which(df_simu$Number_effect ==i  )] )
  
  mean_T1_ps    [h]<- mean(df_simu$t1_ps[which(df_simu$Number_effect ==i )] )
  
  
  mean_cs_size_nps[h] <- tbrm(df_simu$cs_size_nps[which(df_simu$Number_effect ==i )] )
  mean_cs_size_ps[h]  <- tbrm(df_simu$cs_size_ps[which(df_simu$Number_effect ==i )] )
  
  mean_purity_nps[h] <-  mean(df_simu$mean_purity_nps[which(df_simu$Number_effect ==i )] )
  mean_purity_ps[h]  <- mean(df_simu$mean_purity_ps[which(df_simu$Number_effect ==i )] )
  
  
  which_L       [h] <- i
  
  
  h <- h+1
  
  
  
  
}

final_df1 <- rbind(mean_power_nps ,    mean_power_ps,
                   mean_T1_nps ,mean_T1_ps ,
                   mean_cs_size_nps,
                   mean_cs_size_ps,
                   mean_purity_nps,
                   mean_purity_ps,
                   which_L  )

t_names <-rownames(final_df1)
final_df1 <- data.frame(t(final_df1))
colnames(final_df1) <- t_names



df_plot <- data.frame(power=c(final_df1$mean_power_ps, final_df1$mean_power_nps),
                      T1_error =c(final_df1$mean_T1_ps, final_df1$mean_T1_nps),
                      cs_size= c(final_df1$mean_cs_size_ps, final_df1$mean_cs_size_nps),
                      mean_purity =  c(final_df1$mean_purity_ps, final_df1$mean_purity_nps),
                      prior=as.factor(c(rep( "SPSP", nrow(final_df1)),rep( "ISP", nrow(final_df1)))),
                      L=  factor(rep( final_df1$which_L,2)))
df_plot <- df_plot[complete.cases(df_plot),]
sd_error_bin_up <- function(p,n_rep=100){
  c(1.96*sqrt(p*(1-p)/n_rep)+p)
}
sd_error_bin_low <- function(p,n_rep=100){
  c(-1.96*sqrt(p*(1-p)/n_rep)+p)
}

n_rep <- rep(c(table(df_simu$Number_effect)),each=2)

df_plot$pw_er_up <- sd_error_bin_up(df_plot[,1],n_rep)
df_plot$pw_er_low <- sd_error_bin_low(df_plot[,1],n_rep)

df_plot$t1_er_up <- sd_error_bin_up(df_plot[,2],n_rep)
df_plot$t1_er_low <- sd_error_bin_low(df_plot[,2],n_rep)

path <- getwd()
load(paste0(path,
            "/simulation/Simulation_script/script/additional_simualtion_for_fig3_panel_D_susie_pc_calibration_power/check_L_susie_128_block_sd1.RData"))





df_susie= data.frame(do.call( rbind, res))
colnames(df_susie)=c( "Number_effect", "n_cs", "n_effect", "n_false_effect")

df_simu=df_susie
df_simu$power  <- df_simu$n_effect /df_simu$Number_effect 
df_simu$t1  <- 1- df_simu$n_false_effect /(df_simu$n_effect +df_simu$n_false_effect )




mean_power_susie <- rep( NA, max(df_simu$Number_effect) ) 
mean_T1_susie <- rep( NA,  max(df_simu$Number_effect) )

which_L <- rep( NA,  max(df_simu$Number_effect) )

h <- 1
for ( i in unique(df_simu$Number_effect))
{
  
  
  mean_power_susie[h] <- mean(df_simu$power [which(df_simu$Number_effect ==i )] )
  
  mean_T1_susie   [h] <- mean(df_simu$t1 [which(df_simu$Number_effect ==i  )] )
  
  
  which_L       [h] <- i
  
  
  h <- h+1
  
  
  
  
}

mean_T1_susie
mean_power_susie


df= data.frame(power= mean_power_susie,
               T1_error= mean_T1_susie,
               
               L= unique(df_simu$Number_effect))


plot( df$L, df$Power)
plot( df$L, df$T1)


sd_error_bin_up <- function(p,n_rep=100){
  c(1.96*sqrt(p*(1-p)/n_rep)+p)
}
sd_error_bin_low <- function(p,n_rep=100){
  c(-1.96*sqrt(p*(1-p)/n_rep)+p)
}

n_rep <- rep(c(table(df_simu$Number_effect)),each=1)

df$pw_er_up <- sd_error_bin_up(df[,1],n_rep)
df$pw_er_low <- sd_error_bin_low(df[,1],n_rep)

df$t1_er_up <- sd_error_bin_up(df[,2],n_rep)
df$t1_er_low <- sd_error_bin_low(df[,2],n_rep)

df$L=as.factor(df$L)
df$prior=  rep( "SuSiE", nrow(df))


dfp =  rbind ( df_plot[,c("power",
                          "T1_error", 
                          "prior",
                          "L",
                          "pw_er_up",
                          "pw_er_low",
                          "t1_er_up",
                          "t1_er_low") ],
               df     [,c("power",
                          "T1_error", 
                          "prior",
                          "L",
                          "pw_er_up",
                          "pw_er_low",
                          "t1_er_up",
                          "t1_er_low") ])



library(ggplot2)

P21  <- ggplot(dfp, aes( x=L, y= power, col=prior))+
  geom_point(position=position_dodge(.9),size=2)+
  geom_errorbar(aes(ymin=pw_er_low, ymax=pw_er_up), width=.2,
                position=position_dodge(.9) )+
  ylim(c(0 ,1))+
  ylab("Power")+
  theme_linedraw()+
  scale_color_manual(values = colors) +
  scale_x_discrete(labels =c(1,"",3,"",5,"",7,"",9,"",11,"",13,"",15,"",17,"",19,""))

P21 

P21_t1= ggplot(dfp, aes( x=L, y= T1_error, col=prior))+
  
  geom_hline(yintercept = 0.95, color="black", linewidth=2)+
  geom_point(position=position_dodge(.9),size=2)+
  geom_errorbar(aes(ymin=t1_er_low, ymax=t1_er_up), width=.2,
                position=position_dodge(.9) )+
  ylim(c(0.7,1.04))+
  ylab("Coverage")+
  #ggtitle("Coverage PVE=10%") +
  theme_linedraw()+
  scale_color_manual(values = colors)+
  scale_x_discrete(labels =c(1,"",3,"",5,"",7,"",9,"",11,"",13,"",15,"",17,"",19,""))
 
P21_t1

#distance decay -----
path <- getwd()
load(paste( path,"/simulation/Simulation_results/dist_decay_L_accuracy_128_sd1.RData", 
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
#df_simu <- df_simu[which(df_simu$Number_effect %in% c(1,2,4,8,12,16)),]
df_simu <- df_simu[-which(df_simu$cs_size_ps>10),]
library(dplyr)
df_simu$power_nps <- df_simu$n_effect_nps/df_simu$Number_effect
df_simu$power_ps <- df_simu$n_effect_ps/df_simu$Number_effect
df_simu$t1_nps <- 1- df_simu$n_false_effect_nps/(df_simu$n_effect_nps+df_simu$n_false_effect_nps)
df_simu$t1_ps <- 1- df_simu$n_false_effect_ps/(df_simu$n_effect_ps+df_simu$n_false_effect_ps)



mean_power_nps <- rep( NA, max(df_simu$Number_effect) )
mean_power_ps <- rep( NA, max(df_simu$Number_effect) )
mean_T1_nps <- rep( NA,  max(df_simu$Number_effect) )
mean_T1_ps <- rep( NA,  max(df_simu$Number_effect) )

mean_cs_size_nps <- rep( NA,  max(df_simu$Number_effect) )
mean_cs_size_ps <- rep( NA,  max(df_simu$Number_effect) )

mean_purity_nps <- rep( NA,  max(df_simu$Number_effect) )
mean_purity_ps <- rep( NA,  max(df_simu$Number_effect) )
which_L <- rep( NA,  max(df_simu$Number_effect) )

h <- 1
for ( i in unique(df_simu$Number_effect))
{
  
  
  mean_power_nps[h] <- mean(df_simu$power_nps[which(df_simu$Number_effect ==i )] )
  mean_power_ps [h]<- mean(df_simu$power_ps[which(df_simu$Number_effect ==i)  ] )
  
  mean_T1_nps   [h] <- mean(df_simu$t1_nps[which(df_simu$Number_effect ==i  )] )
  
  mean_T1_ps    [h]<- mean(df_simu$t1_ps[which(df_simu$Number_effect ==i )] )
  
  
  mean_cs_size_nps[h] <- tbrm(df_simu$cs_size_nps[which(df_simu$Number_effect ==i )] )
  mean_cs_size_ps[h]  <- tbrm(df_simu$cs_size_ps[which(df_simu$Number_effect ==i )] )
  
  mean_purity_nps[h] <-  mean(df_simu$mean_purity_nps[which(df_simu$Number_effect ==i )] )
  mean_purity_ps[h]  <- mean(df_simu$mean_purity_ps[which(df_simu$Number_effect ==i )] )
  
  
  which_L       [h] <- i
  
  
  h <- h+1
  
  
  
  
}

final_df1 <- rbind(mean_power_nps ,    mean_power_ps,
                   mean_T1_nps ,mean_T1_ps ,
                   mean_cs_size_nps,
                   mean_cs_size_ps,
                   mean_purity_nps,
                   mean_purity_ps,
                   which_L  )

t_names <-rownames(final_df1)
final_df1 <- data.frame(t(final_df1))
colnames(final_df1) <- t_names



df_plot <- data.frame(power=c(final_df1$mean_power_ps, final_df1$mean_power_nps),
                      T1_error =c(final_df1$mean_T1_ps, final_df1$mean_T1_nps),
                      cs_size= c(final_df1$mean_cs_size_ps, final_df1$mean_cs_size_nps),
                      mean_purity =  c(final_df1$mean_purity_ps, final_df1$mean_purity_nps),
                      prior=as.factor(c(rep( "SPSP", nrow(final_df1)),rep( "ISP", nrow(final_df1)))),
                      L=  factor(rep( final_df1$which_L,2)))
df_plot <- df_plot[complete.cases(df_plot),]
sd_error_bin_up <- function(p,n_rep=100){
  c(1.96*sqrt(p*(1-p)/n_rep)+p)
}
sd_error_bin_low <- function(p,n_rep=100){
  c(-1.96*sqrt(p*(1-p)/n_rep)+p)
}

n_rep <- rep(c(table(df_simu$Number_effect)),each=2)

df_plot$pw_er_up <- sd_error_bin_up(df_plot[,1],n_rep)
df_plot$pw_er_low <- sd_error_bin_low(df_plot[,1],n_rep)

df_plot$t1_er_up <- sd_error_bin_up(df_plot[,2],n_rep)
df_plot$t1_er_low <- sd_error_bin_low(df_plot[,2],n_rep)

path <- getwd()
load(paste0(path,
            "/simulation/Simulation_script/script/additional_simualtion_for_fig3_panel_D_susie_pc_calibration_power/check_L_susie_128_distdecay_sd1.RData"))

 
 


df_susie= data.frame(do.call( rbind, res))
colnames(df_susie)=c( "Number_effect", "n_cs", "n_effect", "n_false_effect")

df_simu=df_susie
df_simu$power  <- df_simu$n_effect /df_simu$Number_effect 
df_simu$t1  <- 1- df_simu$n_false_effect /(df_simu$n_effect +df_simu$n_false_effect )
 



mean_power_susie <- rep( NA, max(df_simu$Number_effect) ) 
mean_T1_susie <- rep( NA,  max(df_simu$Number_effect) )
 
which_L <- rep( NA,  max(df_simu$Number_effect) )

h <- 1
for ( i in unique(df_simu$Number_effect))
{
  
  
  mean_power_susie[h] <- mean(df_simu$power [which(df_simu$Number_effect ==i )] )
    
  mean_T1_susie   [h] <- mean(df_simu$t1 [which(df_simu$Number_effect ==i  )] )
  
  
  which_L       [h] <- i
  
  
  h <- h+1
  
  
  
  
}

mean_T1_susie
mean_power_susie


df= data.frame(power= mean_power_susie,
               T1_error= mean_T1_susie,
                
                L= unique(df_simu$Number_effect))


plot( df$L, df$Power)
plot( df$L, df$T1)

 
sd_error_bin_up <- function(p,n_rep=100){
  c(1.96*sqrt(p*(1-p)/n_rep)+p)
}
sd_error_bin_low <- function(p,n_rep=100){
  c(-1.96*sqrt(p*(1-p)/n_rep)+p)
}

n_rep <- rep(c(table(df_simu$Number_effect)),each=1)

df$pw_er_up <- sd_error_bin_up(df[,1],n_rep)
df$pw_er_low <- sd_error_bin_low(df[,1],n_rep)

df$t1_er_up <- sd_error_bin_up(df[,2],n_rep)
df$t1_er_low <- sd_error_bin_low(df[,2],n_rep)

df$L=as.factor(df$L)
df$prior=  rep( "SuSiE", nrow(df))


dfp =  rbind ( df_plot[,c("power",
                          "T1_error", 
                          "prior",
                          "L",
                          "pw_er_up",
                          "pw_er_low",
                          "t1_er_up",
                          "t1_er_low") ],
               df     [,c("power",
                         "T1_error", 
                         "prior",
                         "L",
                         "pw_er_up",
                         "pw_er_low",
                         "t1_er_up",
                         "t1_er_low") ])



library(ggplot2)

P31  <- ggplot(dfp, aes( x=L, y= power, col=prior))+
  geom_point(position=position_dodge(.9),size=2)+
  geom_errorbar(aes(ymin=pw_er_low, ymax=pw_er_up), width=.2,
                position=position_dodge(.9) )+
  ylim(c(0 ,1))+
  ylab("Power")+
  theme_linedraw()+
  scale_color_manual(values = colors) +
  scale_x_discrete(labels =c(1,"",3,"",5,"",7,"",9,"",11,"",13,"",15,"",17,"",19,""))

P31 
P31_t1  <-ggplot(dfp, aes( x=L, y= T1_error, col=prior))+
  
  geom_hline(yintercept = 0.95, color="black", linewidth=2)+
  geom_point(position=position_dodge(.9),size=2)+
  geom_errorbar(aes(ymin=t1_er_low, ymax=t1_er_up), width=.2,
                position=position_dodge(.9) )+
  ylim(c(0.7,1.04))+
  ylab("Coverage")+
  #ggtitle("Coverage PVE=10%") +
  theme_linedraw()+
  scale_color_manual(values = colors)+
  scale_x_discrete(labels =c(1,"",3,"",5,"",7,"",9,"",11,"",13,"",15,"",17,"",19,""))
P31_t1

library( gridExtra)

# Create text labels for the column titles
titles <- list(   textGrob(label = "Gaussian functional", gp = gpar(fontsize = 16, fontface = "bold")),
                  textGrob(label = "WGBS block", gp = gpar(fontsize = 16, fontface = "bold")),
                  textGrob(label = "WGBS distance decay", gp = gpar(fontsize = 16, fontface = "bold"))
)  

legend_plot <- get_legend(P31 ) # Extract all legend components

# Extract the legend from one of the plots
legend_plot <- get_legend(P31)+ theme(legend.position = "bottom",
                                      legend.justification = "center",
                                      
                                      legend.title=element_blank())
grid.arrange(P11,P21,P31, ncol=3)



# Arrange the grid with correct alignment
P_out <- grid.arrange(
  arrangeGrob(
      titles[[1]], titles[[2]],titles[[3]],  # Empty space for alignment
    ncol = 3, widths = c(1, 1, 1)
  ),
  arrangeGrob(
      P11 + theme(legend.position = "none"),
      P21 + theme(legend.position = "none"),
      P31 + theme(legend.position = "none"),
    ncol = 3,widths = c(1, 1, 1)
  ),
 # arrangeGrob(legend_plot),
  heights = c(0.08, 1 )
)

# Display the final plot
P_out


save_path=  paste0(getwd(),
                   "/plot/fig3_separate_panel/"
)
ggsave(P_out , file=paste0(save_path,"power.pdf"),
       width = 29.7,
       height = 10.5,
       units = "cm"
)


# Arrange the grid with correct alignment
P_out <- grid.arrange(
  arrangeGrob(
    titles[[1]], titles[[2]],titles[[3]],  # Empty space for alignment
    ncol = 3, widths = c(1, 1, 1)
  ),
  arrangeGrob(
    P11_t1 + theme(legend.position = "none"),
    P21_t1 + theme(legend.position = "none"),
    P31_t1 + theme(legend.position = "none"),
    ncol = 3,widths = c(1, 1, 1)
  ),
  #arrangeGrob(legend_plot),
  heights = c(0.08, 1 )
)

# Display the final plot
P_out


save_path=  paste0(getwd(),
                   "/plot/fig3_separate_panel/"
)
ggsave(P_out , file=paste0(save_path,"coverage.pdf"),
       width = 29.7,
       height =  10.5,
       units = "cm"
)

