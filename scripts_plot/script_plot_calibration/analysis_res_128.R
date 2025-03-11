library(dplR)
### For PVE=10% ----

path <- getwd()
load(paste( path,"/simulation/Simulation_results/check_L_accuracy_128_sd1.RData", 
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



mean_power_nps <- rep( NA, max(df_simu$Number_effect) )
mean_power_ps <- rep( NA, max(df_simu$Number_effect) )
mean_T1_nps <- rep( NA,  max(df_simu$Number_effect) )
mean_T1_ps <- rep( NA,  max(df_simu$Number_effect) )

mean_cs_size_nps <- rep( NA,  max(df_simu$Number_effect) )
mean_cs_size_ps  <- rep( NA,  max(df_simu$Number_effect) )

var_cs_size_nps <- rep( NA,  max(df_simu$Number_effect) )
var_cs_size_ps  <- rep( NA,  max(df_simu$Number_effect) )


mean_purity_nps <- rep( NA,  max(df_simu$Number_effect) )
mean_purity_ps <- rep( NA,  max(df_simu$Number_effect) )

var_purity_nps <- rep( NA,  max(df_simu$Number_effect) )
var_purity_ps <- rep( NA,  max(df_simu$Number_effect) )

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
  
  
  
  var_cs_size_nps[h] <- var(df_simu$cs_size_nps[which(df_simu$Number_effect ==i )] )
  var_cs_size_ps[h]  <- var(df_simu$cs_size_ps[which(df_simu$Number_effect ==i )] )
  
  var_purity_nps[h] <-  var(df_simu$mean_purity_nps[which(df_simu$Number_effect ==i )] )
  var_purity_ps[h]  <- var(df_simu$mean_purity_ps[which(df_simu$Number_effect ==i )] )
  
  
  which_L       [h] <- i
  
  
  h <- h+1
  
  
  
  
}

final_df1 <- rbind(mean_power_nps ,    mean_power_ps,
                   mean_T1_nps ,mean_T1_ps ,
                   mean_cs_size_nps,
                   mean_cs_size_ps,
                   mean_purity_nps,
                   mean_purity_ps,
                   var_cs_size_nps,
                   var_cs_size_ps,
                   var_purity_nps,
                   var_purity_ps,
                   which_L  )

t_names <-rownames(final_df1)
final_df1 <- data.frame(t(final_df1))
colnames(final_df1) <- t_names



df_plot <- data.frame(power=c(final_df1$mean_power_ps, final_df1$mean_power_nps),
                      T1_error =c(final_df1$mean_T1_ps, final_df1$mean_T1_nps),
                      cs_size= c(final_df1$mean_cs_size_ps, final_df1$mean_cs_size_nps),
                      var_size= c(final_df1$var_cs_size_ps, final_df1$var_cs_size_nps),
                      
                      mean_purity =  c(final_df1$mean_purity_ps, final_df1$mean_purity_nps),
                      var_purity =  c(final_df1$var_purity_ps, final_df1$var_purity_nps),
                      
                      prior=as.factor(c(rep( "SPSP", nrow(final_df1)),rep( "ISP", nrow(final_df1)))),
                      L=  factor(rep( final_df1$which_L,2)))


df_plot <- df_plot[complete.cases(df_plot),]
n_rep <- rep(c(table(df_simu$Number_effect)),each=2)

df_plot$cs_er_up  <-  df_plot$cs_size+ 1.96 *sqrt(df_plot$var_size)/sqrt(100)
df_plot$cs_er_low <-  df_plot$cs_size- 1.96 *sqrt(df_plot$var_size)/sqrt(100)

df_plot$purity_er_up <- df_plot$mean_purity +  1.96 *sqrt(df_plot$var_purity)/sqrt(100)
df_plot$purity_er_low <- df_plot$mean_purity -  1.96 *sqrt(df_plot$var_purity)/sqrt(100)

library(ggplot2)


P1_p <- ggplot(df_plot, aes( x=L, y= power, col=prior))+
  geom_point(position=position_dodge(.9),size=1.2)+
  #geom_errorbar(aes(ymin=pw_er_low, ymax=pw_er_up), width=.2,
  #              position=position_dodge(.9) )+
  ylim(c(0 ,1))+
  ylab("Power")+
  ggtitle("PVE=10%") +
  xlab("Number of effects")+
  scale_color_manual(values=colors)+
  theme_linedraw()
P1_p

P1_t1 <- ggplot(df_plot, aes( x=L, y= T1_error, col=prior))+
  geom_point(position=position_dodge(.9),size=1.2)+
  # geom_errorbar(aes(ymin=t1_er_low, ymax=t1_er_up), width=.2,
  #              position=position_dodge(.9) )+
  geom_hline(yintercept = 0.95)+
  ylim(c(0.8,1.01))+
  ylab("Coverage")+
  xlab("Number of effects")+
  scale_color_manual(values=colors)+
  theme_linedraw()
#ggtitle("Coverage PVE=10%") + 

P1_t1


P1_cs  <- ggplot(df_plot, aes( x=L, y= cs_size, col=prior))+
  geom_point(position=position_dodge(.9),size=1.2)+
  geom_errorbar(aes(ymin=cs_er_low, ymax=cs_er_up), width=.2,
                position=position_dodge(.9) )+
  ylim(c(1,17))+
  ylab("Mean CS size")+
  xlab("Number of effects")+
  scale_color_manual(values=colors)+
  theme_linedraw()
#ggtitle("CS size PVE=10%") + 

P1_cs


P1_pur  <- ggplot(df_plot, aes( x=L, y= mean_purity, col=prior))+
  geom_point(position=position_dodge(.9),size=1.2)+
  geom_errorbar(aes(ymin=purity_er_low, ymax=purity_er_up), width=.2,
                position=position_dodge(.9) )+
  ylim(c(0.9,1))+
  
  ylab("Median purity")+
  xlab("Number of effects")+
  scale_color_manual(values=colors)+
  theme_linedraw()
#ggtitle("Purity size PVE=10%") + 

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
mean_cs_size_ps  <- rep( NA,  max(df_simu$Number_effect) )

var_cs_size_nps <- rep( NA,  max(df_simu$Number_effect) )
var_cs_size_ps  <- rep( NA,  max(df_simu$Number_effect) )


mean_purity_nps <- rep( NA,  max(df_simu$Number_effect) )
mean_purity_ps <- rep( NA,  max(df_simu$Number_effect) )

var_purity_nps <- rep( NA,  max(df_simu$Number_effect) )
var_purity_ps <- rep( NA,  max(df_simu$Number_effect) )

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
  
  
  
  var_cs_size_nps[h] <- var(df_simu$cs_size_nps[which(df_simu$Number_effect ==i )] )
  var_cs_size_ps[h]  <- var(df_simu$cs_size_ps[which(df_simu$Number_effect ==i )] )
  
  var_purity_nps[h] <-  var(df_simu$mean_purity_nps[which(df_simu$Number_effect ==i )] )
  var_purity_ps[h]  <- var(df_simu$mean_purity_ps[which(df_simu$Number_effect ==i )] )
  
  
  which_L       [h] <- i
  
  
  h <- h+1
  
  
  
  
}

final_df1 <- rbind(mean_power_nps ,    mean_power_ps,
                   mean_T1_nps ,mean_T1_ps ,
                   mean_cs_size_nps,
                   mean_cs_size_ps,
                   mean_purity_nps,
                   mean_purity_ps,
                   var_cs_size_nps,
                   var_cs_size_ps,
                   var_purity_nps,
                   var_purity_ps,
                   which_L  )

t_names <-rownames(final_df1)
final_df1 <- data.frame(t(final_df1))
colnames(final_df1) <- t_names



df_plot <- data.frame(power=c(final_df1$mean_power_ps, final_df1$mean_power_nps),
                      T1_error =c(final_df1$mean_T1_ps, final_df1$mean_T1_nps),
                      cs_size= c(final_df1$mean_cs_size_ps, final_df1$mean_cs_size_nps),
                      var_size= c(final_df1$var_cs_size_ps, final_df1$var_cs_size_nps),
                      
                      mean_purity =  c(final_df1$mean_purity_ps, final_df1$mean_purity_nps),
                      var_purity =  c(final_df1$var_purity_ps, final_df1$var_purity_nps),
                      
                      prior=as.factor(c(rep( "SPSP", nrow(final_df1)),rep( "ISP", nrow(final_df1)))),
                      L=  factor(rep( final_df1$which_L,2)))

df_plot <- df_plot[complete.cases(df_plot),]

n_rep <- rep(c(table(df_simu$Number_effect)),each=2)

df_plot$cs_er_up  <-  df_plot$cs_size+ 1.96 *sqrt(df_plot$var_size)/sqrt(100)
df_plot$cs_er_low <-  df_plot$cs_size- 1.96 *sqrt(df_plot$var_size)/sqrt(100)

df_plot$purity_er_up <- df_plot$mean_purity +  1.96 *sqrt(df_plot$var_purity)/sqrt(100)
df_plot$purity_er_low <- df_plot$mean_purity -  1.96 *sqrt(df_plot$var_purity)/sqrt(100)

library(ggplot2)

P2_p <- ggplot(df_plot, aes( x=L, y= power, col=prior))+
  geom_point(position=position_dodge(.9),size=1.2)+
  #geom_errorbar(aes(ymin=pw_er_low, ymax=pw_er_up), width=.2,
  #             position=position_dodge(.9) )+
  ylim(c(0 ,1))+
  ylab(" ")+
  ggtitle("PVE=20%") +
  xlab("Number of effects")+
  scale_color_manual(values=colors)+
  theme_linedraw()
P2_p

P2_t1 <- ggplot(df_plot, aes( x=L, y= T1_error, col=prior))+
  geom_point(position=position_dodge(.9),size=1.2)+
  #geom_errorbar(aes(ymin=t1_er_low, ymax=t1_er_up), width=.2,
  #              position=position_dodge(.9) )+
  geom_hline(yintercept = 0.95)+
  ylim(c(0.8,1.01))+
  ylab(" ")+
  xlab("Number of effects")+
  scale_color_manual(values=colors)+
  theme_linedraw()
#ggtitle("Coverage PVE=20%") + 

P2_t1


P2_cs <- ggplot(df_plot, aes( x=L, y= cs_size, col=prior))+
  geom_point(position=position_dodge(.9),size=1.2)+
  geom_errorbar(aes(ymin=cs_er_low, ymax=cs_er_up), width=.2,
                position=position_dodge(.9) )+
  ylim(c(1,17))+
  ylab(" ")+
  xlab("Number of effects")+
  scale_color_manual(values=colors)+
  theme_linedraw()
#ggtitle("CS size PVE=20%") + 

P2_cs


P2_pur  <- ggplot(df_plot, aes( x=L, y= mean_purity, col=prior))+
  geom_point(position=position_dodge(.9),size=1.2)+
  geom_errorbar(aes(ymin=purity_er_low, ymax=purity_er_up), width=.2,
                position=position_dodge(.9) )+
  ylim(c(0.9,1))+
  
  ylab(" ")+
  xlab("Number of effects")+
  scale_color_manual(values=colors)+
  theme_linedraw()
#ggtitle("Purity size PVE=20%") + 

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
mean_cs_size_ps  <- rep( NA,  max(df_simu$Number_effect) )

var_cs_size_nps <- rep( NA,  max(df_simu$Number_effect) )
var_cs_size_ps  <- rep( NA,  max(df_simu$Number_effect) )


mean_purity_nps <- rep( NA,  max(df_simu$Number_effect) )
mean_purity_ps <- rep( NA,  max(df_simu$Number_effect) )

var_purity_nps <- rep( NA,  max(df_simu$Number_effect) )
var_purity_ps <- rep( NA,  max(df_simu$Number_effect) )

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
  
  
  
  var_cs_size_nps[h] <- var(df_simu$cs_size_nps[which(df_simu$Number_effect ==i )] )
  var_cs_size_ps[h]  <- var(df_simu$cs_size_ps[which(df_simu$Number_effect ==i )] )
  
  var_purity_nps[h] <-  var(df_simu$mean_purity_nps[which(df_simu$Number_effect ==i )] )
  var_purity_ps[h]  <- var(df_simu$mean_purity_ps[which(df_simu$Number_effect ==i )] )
  
  
  which_L       [h] <- i
  
  
  h <- h+1
  
  
  
  
}

final_df1 <- rbind(mean_power_nps ,    mean_power_ps,
                   mean_T1_nps ,mean_T1_ps ,
                   mean_cs_size_nps,
                   mean_cs_size_ps,
                   mean_purity_nps,
                   mean_purity_ps,
                   var_cs_size_nps,
                   var_cs_size_ps,
                   var_purity_nps,
                   var_purity_ps,
                   which_L  )

t_names <-rownames(final_df1)
final_df1 <- data.frame(t(final_df1))
colnames(final_df1) <- t_names



df_plot <- data.frame(power=c(final_df1$mean_power_ps, final_df1$mean_power_nps),
                      T1_error =c(final_df1$mean_T1_ps, final_df1$mean_T1_nps),
                      cs_size= c(final_df1$mean_cs_size_ps, final_df1$mean_cs_size_nps),
                      var_size= c(final_df1$var_cs_size_ps, final_df1$var_cs_size_nps),
                      
                      mean_purity =  c(final_df1$mean_purity_ps, final_df1$mean_purity_nps),
                      var_purity =  c(final_df1$var_purity_ps, final_df1$var_purity_nps),
                      
                      prior=as.factor(c(rep( "SPSP", nrow(final_df1)),rep( "ISP", nrow(final_df1)))),
                      L=  factor(rep( final_df1$which_L,2)))

df_plot <- df_plot[complete.cases(df_plot),]

n_rep <- rep(c(table(df_simu$Number_effect)),each=2)

df_plot$cs_er_up  <-  df_plot$cs_size+ 1.96 *sqrt(df_plot$var_size)/sqrt(100)
df_plot$cs_er_low <-  df_plot$cs_size- 1.96 *sqrt(df_plot$var_size)/sqrt(100)

df_plot$purity_er_up <- df_plot$mean_purity +  1.96 *sqrt(df_plot$var_purity)/sqrt(100)
df_plot$purity_er_low <- df_plot$mean_purity -  1.96 *sqrt(df_plot$var_purity)/sqrt(100)

library(ggplot2)

P3_p <- ggplot(df_plot, aes( x=L, y= power, col=prior))+
  geom_point(position=position_dodge(.9),size=1.2)+
  #geom_errorbar(aes(ymin=pw_er_low, ymax=pw_er_up), width=.2,
  #              position=position_dodge(.9) )+
  ylim(c(0 ,1))+
  ylab(" ")+
  ggtitle("PVE=30%") +
  xlab("Number of effects")+
  scale_color_manual(values=colors)+
  theme_linedraw()
P3_p

P3_t1 <- ggplot(df_plot, aes( x=L, y= T1_error, col=prior))+
  geom_point(position=position_dodge(.9),size=1.2)+
  #geom_errorbar(aes(ymin=t1_er_low, ymax=t1_er_up), width=.2,
  #              position=position_dodge(.9) )+
  geom_hline(yintercept = 0.95)+
  ylim(c(0.8,1.01))+
  ylab(" ")+
  xlab("Number of effects")+
  scale_color_manual(values=colors)+
  theme_linedraw()
#ggtitle("Coverage PVE=30%") + 

P3_t1


P3_cs <- ggplot(df_plot, aes( x=L, y= cs_size, col=prior))+
  geom_point(position=position_dodge(.9),size=1.2)+
  geom_errorbar(aes(ymin=cs_er_low, ymax=cs_er_up), width=.2,
                position=position_dodge(.9) )+
  ylim(c(1,17))+
  ylab(" ")+
  xlab("Number of effects")+
  scale_color_manual(values=colors)+
  theme_linedraw()
#ggtitle("PVE=30%") +


P3_cs


P3_pur  <- ggplot(df_plot, aes( x=L, y= mean_purity, col=prior))+
  geom_point(position=position_dodge(.9),size=1.2)+
  geom_errorbar(aes(ymin=purity_er_low, ymax=purity_er_up), width=.2,
                position=position_dodge(.9) )+
  ylim(c(0.9,1))+
  
  ylab(" ")+
  xlab("Number of effects")+
  scale_color_manual(values=colors)+
  theme_linedraw()
#ggtitle("Purity size PVE=30%") + 

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
mean_cs_size_ps  <- rep( NA,  max(df_simu$Number_effect) )

var_cs_size_nps <- rep( NA,  max(df_simu$Number_effect) )
var_cs_size_ps  <- rep( NA,  max(df_simu$Number_effect) )


mean_purity_nps <- rep( NA,  max(df_simu$Number_effect) )
mean_purity_ps <- rep( NA,  max(df_simu$Number_effect) )

var_purity_nps <- rep( NA,  max(df_simu$Number_effect) )
var_purity_ps <- rep( NA,  max(df_simu$Number_effect) )

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
  
  
  
  var_cs_size_nps[h] <- var(df_simu$cs_size_nps[which(df_simu$Number_effect ==i )] )
  var_cs_size_ps[h]  <- var(df_simu$cs_size_ps[which(df_simu$Number_effect ==i )] )
  
  var_purity_nps[h] <-  var(df_simu$mean_purity_nps[which(df_simu$Number_effect ==i )] )
  var_purity_ps[h]  <- var(df_simu$mean_purity_ps[which(df_simu$Number_effect ==i )] )
  
  
  which_L       [h] <- i
  
  
  h <- h+1
  
  
  
  
}

final_df1 <- rbind(mean_power_nps ,    mean_power_ps,
                   mean_T1_nps ,mean_T1_ps ,
                   mean_cs_size_nps,
                   mean_cs_size_ps,
                   mean_purity_nps,
                   mean_purity_ps,
                   var_cs_size_nps,
                   var_cs_size_ps,
                   var_purity_nps,
                   var_purity_ps,
                   which_L  )

t_names <-rownames(final_df1)
final_df1 <- data.frame(t(final_df1))
colnames(final_df1) <- t_names



df_plot <- data.frame(power=c(final_df1$mean_power_ps, final_df1$mean_power_nps),
                      T1_error =c(final_df1$mean_T1_ps, final_df1$mean_T1_nps),
                      cs_size= c(final_df1$mean_cs_size_ps, final_df1$mean_cs_size_nps),
                      var_size= c(final_df1$var_cs_size_ps, final_df1$var_cs_size_nps),
                      
                      mean_purity =  c(final_df1$mean_purity_ps, final_df1$mean_purity_nps),
                      var_purity =  c(final_df1$var_purity_ps, final_df1$var_purity_nps),
                      
                      prior=as.factor(c(rep( "SPSP", nrow(final_df1)),rep( "ISP", nrow(final_df1)))),
                      L=  factor(rep( final_df1$which_L,2)))

df_plot <- df_plot[complete.cases(df_plot),]

n_rep <- rep(c(table(df_simu$Number_effect)),each=2)

df_plot$cs_er_up  <-  df_plot$cs_size+ 1.96 *sqrt(df_plot$var_size)/sqrt(100)
df_plot$cs_er_low <-  df_plot$cs_size- 1.96 *sqrt(df_plot$var_size)/sqrt(100)

df_plot$purity_er_up <- df_plot$mean_purity +  1.96 *sqrt(df_plot$var_purity)/sqrt(100)
df_plot$purity_er_low <- df_plot$mean_purity -  1.96 *sqrt(df_plot$var_purity)/sqrt(100)

library(ggplot2)

P4_p <- ggplot(df_plot, aes( x=L, y= power, col=prior))+
  geom_point(position=position_dodge(.9),size=1.2)+
  #geom_errorbar(aes(ymin=pw_er_low, ymax=pw_er_up), width=.2,
  #              position=position_dodge(.9) )+
  ylim(c(0 ,1))+
  ylab(" ")+
  ggtitle(" PVE=40%")+
  xlab("Number of effects")+
  scale_color_manual(values=colors)+
  theme_linedraw()
P4_p

P4_t1 <- ggplot(df_plot, aes( x=L, y= T1_error, col=prior))+
  geom_point(position=position_dodge(.9),size=1.2)+
  #geom_errorbar(aes(ymin=t1_er_low, ymax=t1_er_up), width=.2,
  #              position=position_dodge(.9) )+
  geom_hline(yintercept = 0.95)+
  ylim(c(0.8,1.01))+
  ylab(" ")+
  xlab("Number of effects")+
  scale_color_manual(values=colors)+
  theme_linedraw()
#ggtitle("Coverage PVE=40%") + 

P4_t1


P4_cs <- ggplot(df_plot, aes( x=L, y= cs_size, col=prior))+
  geom_point(position=position_dodge(.9),size=1.2)+
  geom_errorbar(aes(ymin=cs_er_low, ymax=cs_er_up), width=.2,
                position=position_dodge(.9) )+
  ylim(c(1,17))+
  ylab(" ")+
  xlab("Number of effects")+
  scale_color_manual(values=colors)+
  theme_linedraw()
#ggtitle("CS size PVE=40%") +


P4_cs


P4_pur  <- ggplot(df_plot, aes( x=L, y= mean_purity, col=prior))+
  geom_point(position=position_dodge(.9),size=1.2)+
  geom_errorbar(aes(ymin=purity_er_low, ymax=purity_er_up), width=.2,
                position=position_dodge(.9) )+
  ylim(c(0.9,1))+
  ylab(" ") +
  xlab("Number of effects")+
  scale_color_manual(values=colors)+
  theme_linedraw()
#ggtitle("Purity size PVE=40%") + 

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
path= getwd( )
# Open a PDF device with A4 landscape dimensions

pdf(paste0 (path, "/plot/simu_128.pdf") , width = 11.69, height = 8.27)  
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

dev.off()
