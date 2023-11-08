library(dplR)
### For PVE=10% ----

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
df_simu <- df_simu[which(df_simu$Number_effect %in% c(1,2,4,8,12,16)),]
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
  c(sqrt(p*(1-p)/n_rep)+p)
}
sd_error_bin_low <- function(p,n_rep=100){
  c(-sqrt(p*(1-p)/n_rep)+p)
}

n_rep <- rep(c(table(df_simu$Number_effect)),each=2)

df_plot$pw_er_up <- sd_error_bin_up(df_plot[,1],n_rep)
df_plot$pw_er_low <- sd_error_bin_low(df_plot[,1],n_rep)

df_plot$t1_er_up <- sd_error_bin_up(df_plot[,2],n_rep)
df_plot$t1_er_low <- sd_error_bin_low(df_plot[,2],n_rep)

library(ggplot2)

P1_p <- ggplot(df_plot, aes( x=L, y= power, col=prior))+
  geom_point(position=position_dodge(.9),size=2)+
  geom_errorbar(aes(ymin=pw_er_low, ymax=pw_er_up), width=.2,
                position=position_dodge(.9) )+
  ylim(c(0 ,1))+
  ylab("Power")+
  ggtitle("PVE=10%") +
  theme_linedraw()
P1_p

P1_t1 <- ggplot(df_plot, aes( x=L, y= T1_error, col=prior))+
  geom_point(position=position_dodge(.9),size=2)+
  geom_errorbar(aes(ymin=t1_er_low, ymax=t1_er_up), width=.2,
                position=position_dodge(.9) )+
  geom_hline(yintercept = 0.95)+
  ylim(c(0.8,1.01))+
  ylab("Coverage")+
  #ggtitle("Coverage PVE=10%") +
  theme_linedraw()

P1_t1


P1_cs  <- ggplot(df_plot, aes( x=L, y= cs_size, col=prior))+
  geom_point(position=position_dodge(.9),size=2)+

  ylim(c(1,17))+
  ylab("CS size")+
  #ggtitle("CS size PVE=10%") +
  theme_linedraw()

P1_cs


P1_pur  <- ggplot(df_plot, aes( x=L, y= mean_purity, col=prior))+
  geom_point(position=position_dodge(.9),size=2)+

  ylim(c(0.9,1))+

  ylab("Purity")+
  #ggtitle("Purity size PVE=10%") +
  theme_linedraw()

P1_pur
### For PVE=20% ------ 
path <- getwd()
load(paste( path,"/simulation/Simulation_results/check_L_accuracy_128_sd2.RData", 
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
  c(sqrt(p*(1-p)/n_rep)+p)
}
sd_error_bin_low <- function(p,n_rep=100){
  c(-sqrt(p*(1-p)/n_rep)+p)
}

n_rep <- rep(c(table(df_simu$Number_effect)),each=2)

df_plot$pw_er_up <- sd_error_bin_up(df_plot[,1],n_rep)
df_plot$pw_er_low <- sd_error_bin_low(df_plot[,1],n_rep)

df_plot$t1_er_up <- sd_error_bin_up(df_plot[,2],n_rep)
df_plot$t1_er_low <- sd_error_bin_low(df_plot[,2],n_rep)

library(ggplot2)

P2_p <- ggplot(df_plot, aes( x=L, y= power, col=prior))+
  geom_point(position=position_dodge(.9),size=2)+
  geom_errorbar(aes(ymin=pw_er_low, ymax=pw_er_up), width=.2,
                position=position_dodge(.9) )+
  ylim(c(0 ,1))+
  ylab(" ")+
  ggtitle("PVE=20%") +
  theme_linedraw()
P2_p

P2_t1 <- ggplot(df_plot, aes( x=L, y= T1_error, col=prior))+
  geom_point(position=position_dodge(.9),size=2)+
  geom_errorbar(aes(ymin=t1_er_low, ymax=t1_er_up), width=.2,
                position=position_dodge(.9) )+
  geom_hline(yintercept = 0.95)+
  ylim(c(0.8,1.01))+
   ylab(" ")+
  #ggtitle("Coverage PVE=20%") +
  theme_linedraw()

P2_t1


P2_cs <- ggplot(df_plot, aes( x=L, y= cs_size, col=prior))+
  geom_point(position=position_dodge(.9),size=2)+

  ylim(c(1,17))+
   ylab(" ")+
  #ggtitle("CS size PVE=20%") +
  theme_linedraw()

P2_cs


P2_pur  <- ggplot(df_plot, aes( x=L, y= mean_purity, col=prior))+
  geom_point(position=position_dodge(.9),size=2)+

  ylim(c(0.9,1))+

  ylab(" ")+
  #ggtitle("Purity size PVE=20%") +
  theme_linedraw()

P2_pur
### For PVE=30% -----

path <- getwd()
load(paste( path,"/simulation/Simulation_results/check_L_accuracy_128_sd3.RData", 
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
  c(sqrt(p*(1-p)/n_rep)+p)
}
sd_error_bin_low <- function(p,n_rep=100){
  c(-sqrt(p*(1-p)/n_rep)+p)
}

n_rep <- rep(c(table(df_simu$Number_effect)),each=2)

df_plot$pw_er_up <- sd_error_bin_up(df_plot[,1],n_rep)
df_plot$pw_er_low <- sd_error_bin_low(df_plot[,1],n_rep)

df_plot$t1_er_up <- sd_error_bin_up(df_plot[,2],n_rep)
df_plot$t1_er_low <- sd_error_bin_low(df_plot[,2],n_rep)
library(ggplot2)

P3_p <- ggplot(df_plot, aes( x=L, y= power, col=prior))+
  geom_point(position=position_dodge(.9),size=2)+
  geom_errorbar(aes(ymin=pw_er_low, ymax=pw_er_up), width=.2,
                position=position_dodge(.9) )+
  ylim(c(0 ,1))+
   ylab(" ")+
   ggtitle("PVE=30%")  +
  theme_linedraw()
P3_p

P3_t1 <- ggplot(df_plot, aes( x=L, y= T1_error, col=prior))+
  geom_point(position=position_dodge(.9),size=2)+
  geom_errorbar(aes(ymin=t1_er_low, ymax=t1_er_up), width=.2,
                position=position_dodge(.9) )+
  geom_hline(yintercept = 0.95)+
  ylim(c(0.8,1.01))+
   ylab(" ")+
  #ggtitle("Coverage PVE=30%") +
  theme_linedraw()

P3_t1


P3_cs <- ggplot(df_plot, aes( x=L, y= cs_size, col=prior))+
  geom_point(position=position_dodge(.9),size=2)+

  ylim(c(1,17))+
   ylab(" ")+
  #ggtitle("PVE=30%") +
  theme_linedraw()

P3_cs


P3_pur  <- ggplot(df_plot, aes( x=L, y= mean_purity, col=prior))+
  geom_point(position=position_dodge(.9),size=2)+

  ylim(c(0.9,1))+

   ylab(" ")+
  #ggtitle("Purity size PVE=30%") +
  theme_linedraw()

P3_pur
### For PVE=40% ------ 
path <- getwd()
load(paste( path,"/simulation/Simulation_results/check_L_accuracy_128_sd4.RData", 
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

  mean_purity_nps[h] <-  median(df_simu$mean_purity_nps[which(df_simu$Number_effect ==i )] )
  mean_purity_ps[h]  <- median(df_simu$mean_purity_ps[which(df_simu$Number_effect ==i )] )


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
  c(sqrt(p*(1-p)/n_rep)+p)
}
sd_error_bin_low <- function(p,n_rep=100){
  c(-sqrt(p*(1-p)/n_rep)+p)
}

n_rep <- rep(c(table(df_simu$Number_effect)),each=2)

df_plot$pw_er_up <- sd_error_bin_up(df_plot[,1],n_rep)
df_plot$pw_er_low <- sd_error_bin_low(df_plot[,1],n_rep)

df_plot$t1_er_up <- sd_error_bin_up(df_plot[,2],n_rep)
df_plot$t1_er_low <- sd_error_bin_low(df_plot[,2],n_rep)
library(ggplot2)

P4_p <- ggplot(df_plot, aes( x=L, y= power, col=prior))+
  geom_point(position=position_dodge(.9),size=2)+
  geom_errorbar(aes(ymin=pw_er_low, ymax=pw_er_up), width=.2,
                position=position_dodge(.9) )+
  ylim(c(0 ,1))+
   ylab(" ")+
  ggtitle(" PVE=40%") +
  theme_linedraw()
P4_p

P4_t1 <- ggplot(df_plot, aes( x=L, y= T1_error, col=prior))+
  geom_point(position=position_dodge(.9),size=2)+
  geom_errorbar(aes(ymin=t1_er_low, ymax=t1_er_up), width=.2,
                position=position_dodge(.9) )+
  geom_hline(yintercept = 0.95)+
  ylim(c(0.8,1.01))+
   ylab(" ")+
  #ggtitle("Coverage PVE=40%") +
  theme_linedraw()

P4_t1


P4_cs <- ggplot(df_plot, aes( x=L, y= cs_size, col=prior))+
  geom_point(position=position_dodge(.9),size=2)+

  ylim(c(1,17))+
   ylab(" ")+
  #ggtitle("CS size PVE=40%") +
  theme_linedraw()

P4_cs


P4_pur  <- ggplot(df_plot, aes( x=L, y= mean_purity, col=prior))+
  geom_point(position=position_dodge(.9),size=2)+

  ylim(c(0.9,1))+
   ylab(" ")+
  #ggtitle("Purity size PVE=40%") +
  theme_linedraw()

P4_pur

#### gatrhering results

library(gridExtra)

grid.arrange(P1_p,P2_p,P3_p, P4_p, ncol=2)

grid.arrange(P1_t1,P2_t1,P3_t1, P4_t1, ncol=2)

grid.arrange(P1_cs,P2_cs,P3_cs, P4_cs, ncol=2)
grid.arrange(P1_pur,P2_pur,P3_pur, P4_pur, ncol=2)


grid.arrange(P1_p,P2_p,P3_p, P4_p ,
             P1_t1,P2_t1,P3_t1, P4_t1,
             P1_cs,P2_cs,P3_cs, P4_cs,  P1_pur,P2_pur,P3_pur, P4_pur, ncol=4)



library(grid)

library(ggpubr)
ggarrange(P1_p,P2_p,P3_p,P4_p,
          P1_t1,P2_t1,P3_t1,P4_t1,
          P1_cs,P2_cs,P3_cs,P4_cs,
          P1_pur,P2_pur,P3_pur,P4_pur ,
nrow=4,
          ncol=4,common.legend = TRUE, legend="bottom")
x = c(0.01  , 0.5, 0.99,       0.5 , 0.5 , 0.5,0.25 , 0.25 ,0.25,0.75 , 0.75 ,0.75,    0.01  , 0.5, 0.99,    0.01  , 0.5, 0.99  )
y = c(0.525, 0.525, 0.525, 0.05   , 0.5 , 0.99 , 0.05   , 0.5 , 0.99, 0.05   , 0.5 , 0.99, 0.76, 0.76, 0.76, 0.285, 0.285, 0.285)
id = c(1, 1    ,1               ,2   ,2   ,2 ,3         ,3   ,3       ,4   ,4   ,4      ,5   ,5   ,5     ,6   ,6   ,6)


grid.polygon(x,y,id)


