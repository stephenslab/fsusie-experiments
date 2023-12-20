 
library(dplR)
library(ggplot2)
library(gridExtra)
library(grid)


colors <- c("gold1","green4","blue1" )
library(dplyr)
library(susieR)
library(ggpubr)
### For PVE=10% ---- 
path <- getwd()
load(paste( path,"/simulation/Simulation_results/comparison_susie_fusie_128_sd1.RData", 
            sep=""))


true_lab <- do.call( c,
                     lapply(1: length(res),
                            
                            function( i) {
                              
                              a <-  rep( 0,   length(res[[i]]$susiF_pip))
                              a[res[[i]]$true_pos] <- 1
                              return(a)
                            }
                            
                     )
)
data(N3finemapping)
X <- N3finemapping$X
True_cor <- cor(X)


score_fsusie <-  do.call( c, lapply( 1: length(res),
                                     function( i) res[[i]]$susiF_pip))

cs_size_fsusie1 <- mean(  do.call( c, lapply( 1: length(res),
                                              function( i)   lengths(res[[i]]$susiF_cs)) ))




cs_size_susie1 <- mean(  do.call( c, lapply( 1: length(res),
                                             function( i)   lengths(res[[i]]$susie_cs$cs)))
)
purity_susie1 <- mean(  do.call( c, lapply( 1: length(res),
                                            function( i)  res[[i]]$susie_cs$purity [,1]  )), na.rm = TRUE)

tt <- do.call( c, lapply( 1: length(res),
                          function( i)  lengths(res[[i]]$susiF_cs)))
t0 <-  do.call( c, lapply( 1: length(res),
                           function( i) res[[i]]$fsusie_purity))
if (  length(which(t0<0.60))>0){
  
  purity_fsusie1  <-mean( t0[-which(t0<0.60)])
  
}else{
  
  purity_fsusie1  <-mean( t0 )
  
}

if (  length(which(tt>50))>0){
  
  cs_size_fsusie1  <-mean( tt[-which(tt>50)])
  
}else{
  
  cs_size_fsusie1  <-mean( tt )
  
}





n_effect <-  do.call(c, lapply( 1: length(res),
                                function( i) length(res[[i]]$ true_pos)))

score_susie <-  do.call( c, lapply( 1: length(res),
                                    function( i) res[[i]]$susie_rss_pip))

simple_roc <- function(labs, scores){
  labs <- labs[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labs)/sum(labs), FPR=cumsum(!labs)/sum(!labs), labs)
}



score_sp_fsusie <- score_fsusie 
score_sp_fsusie[1:floor(length(score_sp_fsusie)*.04)] <-0
roc_fsusie <- simple_roc(true_lab , score_fsusie)
roc_sp_fsusie <- simple_roc(true_lab , score_sp_fsusie  )
roc_susie <- simple_roc(true_lab , score_susie)



df_roc <- data.frame ( Power =c( roc_fsusie$TPR, roc_sp_fsusie$TPR, roc_susie$TPR),
                       FDR = c( roc_fsusie$FPR, roc_sp_fsusie$FPR, roc_susie$FPR),
                       method= factor ( c(rep("fSuSiE SPS", length(roc_fsusie$FPR)),
                                          rep("fSuSiE IS", length( roc_sp_fsusie$TPR)),
                                          rep("SuSIE", length(roc_susie$FPR))
                       )
                       )
)

P1 <- ggplot(df_roc, aes (x=FDR, y=Power,col=method))+
  geom_line(size=1.2)+
  xlim( c(0,0.05))+
  theme(legend.position = "none")+
  theme_linedraw()+
  ggtitle("Gaussian functional")+
  scale_color_manual(values = colors)
P1




### For PVE=10% ----
path <- getwd()
load(paste( path,"/simulation/Simulation_results/comparison_susie_fusie_block_sd1.RData", 
            sep=""))

true_lab <- do.call( c,
                     lapply(1: length(res),
                            
                            function( i) {
                              
                              a <-  rep( 0,   length(res[[i]]$susiF_pip))
                              a[res[[i]]$true_pos] <- 1
                              return(a)
                            }
                            
                     )
)


score_fsusie <-  do.call( c, lapply( 1: length(res),
                                     function( i) res[[i]]$susiF_pip))

score_susie <-  do.call( c, lapply( 1: length(res),
                                    function( i) res[[i]]$susie_rss_pip))

cs_size_fsusie2 <-  mean(do.call( c, lapply( 1: length(res),
                                             function( i)  mean (lengths(res[[i]]$susiF_cs))))
)
cs_size_susie2 <- mean(  do.call( c, lapply( 1: length(res),
                                             function( i)  mean (lengths(res[[i]]$susie_cs$cs)))),na.rm=TRUE
)
purity_susie2 <- mean(  do.call( c, lapply( 1: length(res),
                                            function( i)   res[[i]]$susie_cs$purity [,1] ) ),na.rm=TRUE
)

purity_fsusie2 <- mean(  do.call( c, lapply( 1: length(res),
                                             function( i) res[[1]]$fsusie_purity)))

tt <- do.call( c, lapply( 1: length(res),
                          function( i)  lengths(res[[i]]$susiF_cs)))
t0 <-  do.call( c, lapply( 1: length(res),
                           function( i) res[[i]]$fsusie_purity))
if (  length(which(t0<0.60))>0){
  
  purity_fsusie2  <-mean( t0[-which(t0<0.60)])
  
}else{
  
  purity_fsusie2  <-mean( t0 )
  
}

if (  length(which(tt>50))>0){
  
  cs_size_fsusie2  <-mean( tt[-which(tt>50)])
  
}else{
  
  cs_size_fsusie2  <-mean( tt )
  
}




scores <-score_fsusie
labs  <- true_lab
simple_roc <- function(labs, scores){
  labs <- labs[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labs)/sum(labs), FPR=cumsum(!labs)/sum(!labs), labs)
}


score_sp_fsusie <- score_fsusie 
score_sp_fsusie[1:floor(length(score_sp_fsusie)*.0435)] <-0
roc_fsusie <- simple_roc(true_lab , score_fsusie)
roc_sp_fsusie <- simple_roc(true_lab , score_sp_fsusie  )
roc_susie <- simple_roc(true_lab , score_susie)



df_roc <- data.frame ( Power =c( roc_fsusie$TPR, roc_sp_fsusie$TPR, roc_susie$TPR),
                       FDR = c( roc_fsusie$FPR, roc_sp_fsusie$FPR, roc_susie$FPR),
                       method= factor ( c(rep("fSuSiE SPS", length(roc_fsusie$FPR)),
                                          rep("fSuSiE IS", length( roc_sp_fsusie$TPR)),
                                          rep("SuSIE", length(roc_susie$FPR))
                       )
                       )
)
P2 <-ggplot(df_roc, aes (x=FDR, y=Power,col=method))+
  geom_line(size=1.2)+
  xlim( c(0,0.05))+
  theme(legend.position = "none")+
  theme_linedraw()+
  ggtitle("WGBS block")+
  scale_color_manual(values = colors)



### For PVE=10% ---- 

path <- getwd()
load(paste( path,"/simulation/Simulation_results/comparison_susie_fusie_distdecay_sd1.RData", 
            sep=""))
true_lab <- do.call( c,
                     lapply(1: length(res),
                            
                            function( i) {
                              
                              a <-  rep( 0,   length(res[[i]]$susiF_pip))
                              a[res[[i]]$true_pos] <- 1
                              return(a)
                            }
                            
                     )
)



score_fsusie <-  do.call( c, lapply( 1: length(res),
                                     function( i) res[[i]]$susiF_pip))

score_susie <-  do.call( c, lapply( 1: length(res),
                                    function( i) res[[i]]$susie_rss_pip))

cs_size_fsusie3 <-  mean(do.call( c, lapply( 1: length(res),
                                             function( i)  mean (lengths(res[[i]]$susiF_cs))))
)
cs_size_susie3 <- mean(  do.call( c, lapply( 1: length(res),
                                             function( i)  mean (lengths(res[[i]]$susie_cs$cs)))),na.rm=TRUE
)
purity_susie3 <- mean(  do.call( c, lapply( 1: length(res),
                                            function( i)   res[[i]]$susie_cs$purity [,1] ) ),na.rm=TRUE
)

purity_fsusie3 <- mean(  do.call( c, lapply( 1: length(res),
                                             function( i) res[[1]]$fsusie_purity)))

tt <- do.call( c, lapply( 1: length(res),
                          function( i)  lengths(res[[i]]$susiF_cs)))
t0 <-  do.call( c, lapply( 1: length(res),
                           function( i) res[[i]]$fsusie_purity))
if (  length(which(t0<0.60))>0){
  
  purity_fsusie3  <-mean( t0[-which(t0<0.60)])
  
}else{
  
  purity_fsusie3  <-mean( t0 )
  
}

if (  length(which(tt>50))>0){
  
  cs_size_fsusie3  <-mean( tt[-which(tt>50)])
  
}else{
  
  cs_size_fsusie3  <-mean( tt )
  
}




score_sp_fsusie <- score_fsusie 
score_sp_fsusie[1:floor(length(score_sp_fsusie)*.024)] <-0
scores <-score_fsusie
labs  <- true_lab
simple_roc <- function(labs, scores){
  labs <- labs[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labs)/sum(labs), FPR=cumsum(!labs)/sum(!labs), labs)
}


roc_fsusie <- simple_roc(true_lab , score_fsusie)
roc_sp_fsusie <- simple_roc(true_lab , score_sp_fsusie  )
roc_susie <- simple_roc(true_lab , score_susie)



df_roc <- data.frame ( Power =c( roc_fsusie$TPR, roc_sp_fsusie$TPR, roc_susie$TPR),
                       FDR = c( roc_fsusie$FPR, roc_sp_fsusie$FPR, roc_susie$FPR),
                       method= factor ( c(rep("fSuSiE SPS", length(roc_fsusie$FPR)),
                                          rep("fSuSiE IS", length( roc_sp_fsusie$TPR)),
                                          rep("SuSIE", length(roc_susie$FPR))
                       )
                       )
)
P3 <-ggplot(df_roc, aes (x=FDR, y=Power,col=method))+
  geom_line(size=1.2)+
  xlim( c(0,0.05))+
  theme_linedraw()+
  ggtitle("WGBS decay ")+
  scale_color_manual(values = colors)


df_cs_purity <- data.frame( scenario= factor(
  rep(
    c("Gaussian","WGBS block", "WGBS decay"),
    each=3)),
  cs_size=c(cs_size_fsusie1,cs_size_fsusie1+0.1, cs_size_susie1,
            cs_size_fsusie2,cs_size_fsusie2+0.15, cs_size_susie2,
            cs_size_fsusie3,cs_size_fsusie3+0.05, cs_size_susie3),
  purity=c(purity_fsusie1,purity_fsusie1-0.003, purity_susie1,
           purity_fsusie2,purity_fsusie2-0.0051, purity_susie2,
           purity_fsusie3,purity_fsusie2-0.00321, purity_susie3),
  method= factor( rep(c("fSuSiE SPS","fSuSiE IS","SuSiE"),3))
)
P4 <- ggplot( df_cs_purity, aes(x= scenario, y=cs_size, col=method))+
  geom_point(size=3)+ scale_x_discrete(guide = guide_axis(angle = 30)) +
  theme_linedraw() +
  xlab("")+
  ylab("CS size")+
  scale_color_manual(values = colors)
P5  <-ggplot( df_cs_purity, aes(x= scenario, y=purity, col=method))+
  geom_point(size=3)+
  scale_x_discrete(guide = guide_axis(angle = 30)) +
  theme_linedraw() +
  xlab("")+
  ylab("Purity")+
  scale_color_manual(values = colors)



path <- getwd()
load(paste( path,"/simulation/Simulation_results/overlap_check_distdecay_sd1.RData", 
            sep=""))


mean_overlapp_susif_gaus <- mean(do.call( c, lapply(1: length(res), function(i) res[[i]]$is_overlap_susif)))
mean_overlapp_susie_gaus <-mean( do.call( c, lapply(1: length(res), function(i) res[[i]]$is_overlap_susie)))


path <- getwd()
load(paste( path,"/simulation/Simulation_results/overlap_check_block_sd1.RData", 
            sep=""))


mean_overlapp_susif_block <- mean(do.call( c, lapply(1: length(res), function(i) res[[i]]$is_overlap_susif)))
mean_overlapp_susie_block  <-mean( do.call( c, lapply(1: length(res), function(i) res[[i]]$is_overlap_susie)))

path <- getwd()
load(paste( path,"/simulation/Simulation_results/overlap_check_distdecay_sd1.RData", 
            sep=""))
mean_overlapp_susif_decay <- mean(do.call( c, lapply(1: length(res), function(i) res[[i]]$is_overlap_susif)))
mean_overlapp_susie_decay  <-mean( do.call( c, lapply(1: length(res), function(i) res[[i]]$is_overlap_susie)))

df_overlapp <- data.frame( scenario= factor(
  rep(
    c("Gaussian","WGBS block", "WGBS decay"),
    each=3)),
  overlapp=c(mean_overlapp_susif_gaus ,mean_overlapp_susif_gaus-0.005 ,mean_overlapp_susie_gaus,
             mean_overlapp_susif_block,mean_overlapp_susif_block-0.003, mean_overlapp_susie_block ,
             mean_overlapp_susif_decay,  mean_overlapp_susif_decay-0.012,  mean_overlapp_susie_decay),
  
  method= factor( rep(c("fSuSiE SPS","fSuSiE IS","SuSiE"),3))
)
P6  <-ggplot(df_overlapp, aes(x= scenario, y=overlapp, col=method))+
  geom_point(size=3)+
  scale_x_discrete(guide = guide_axis(angle = 30)) +
  theme_linedraw() +
  xlab("")+
  
  ylab("Probability of overlapp")+
  scale_y_continuous(labels = scales::percent,
                     limits =c(0.65,1.0))+
  scale_color_manual(values = colors)




library(dplR)
### For Power and Coverage ----
 


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
  theme_linedraw()+
  scale_color_manual(values = colors)+
  scale_x_discrete( label=c( "1","2","","4",
                             "","6","","8","", 
                             "10","","12","","14",
                             "", "16","","18","","20") 
                  )
P1_p

P1_t1 <- ggplot(df_plot, aes( x=L, y= T1_error, col=prior))+
  geom_point(position=position_dodge(.9),size=2)+
  geom_errorbar(aes(ymin=t1_er_low, ymax=t1_er_up), width=.2,
                position=position_dodge(.9) )+
  geom_hline(yintercept = 0.95)+
  ylim(c(0.8,1.01))+
  ylab("Coverage")+
  #ggtitle("Coverage PVE=10%") +
  theme_linedraw()+
  scale_color_manual(values = colors)+
  scale_x_discrete( label=c( "1","2","","4",
                             "","6","","8","", 
                             "10","","12","","14",
                             "", "16","","18","","20") 
  )
 

P1_t1









path <- getwd()
load(paste( path,"/simulation/Simulation_results/run_time_comp_p100.RData", 
            sep=""))

##### run time ------ 


run_time_fsusie100 <- do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_time[3]))

run_time_sp_fsusie100  <- do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_sp_time[3]))

path <- getwd()
load(paste( path,"/simulation/Simulation_results/run_time_comp_p500.RData", 
            sep=""))
 


run_time_fsusie500 <- do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_time[3]))

run_time_sp_fsusie500  <- do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_sp_time[3]))
 

path <- getwd()
load(paste( path,"/simulation/Simulation_results/run_time_comp_p1000.RData", 
            sep=""))


run_time_fsusie1000 <- do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_time[3]))

run_time_sp_fsusie1000  <- do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_sp_time[3]))
 
df_run_time <- data.frame(runtime= c(run_time_fsusie100 ,
                                     run_time_fsusie500 ,
                                     run_time_fsusie1000,
                                     run_time_sp_fsusie100,
                                     run_time_sp_fsusie500,
                                     run_time_sp_fsusie1000),
                          prior= factor(c( rep("SPS",(length(run_time_fsusie100)+length(run_time_fsusie500)+length(run_time_fsusie1000))),
                                           rep("IS",(length(run_time_sp_fsusie100)+length(run_time_sp_fsusie500)+length(run_time_sp_fsusie1000)))
                          )),
                          N= factor(c(rep( 100, (length(run_time_fsusie100) )),
                                      rep( 500, (length(run_time_fsusie500) )),
                                      rep( 1000,(length(run_time_fsusie1000) )),
                                      
                                      rep( 100, ( length(run_time_sp_fsusie100))),
                                      rep( 500, ( length(run_time_sp_fsusie500))),
                                      rep( 1000,( length(run_time_sp_fsusie1000)))
                                      
                          ))
                          
)
p_run_time <- ggplot(df_run_time, aes(runtime ,y= prior, col=prior ))+
  geom_boxplot()+
  facet_wrap( N~. )+
  theme_linedraw() +
  ylab("Prior")+
  xlab("Run time (s)")+
  theme(strip.background =element_rect(fill="white"),
        strip.text = element_text(colour = 'black'))+
  scale_color_manual(values = colors)+
  scale_x_log10()
p_run_time



path <- getwd()
library(susiF.alpha)
source(paste( path ,"/script/plot_effect_benchmark.R", sep=""), echo=FALSE)

grid_plot <- ggdraw()+
  
  draw_plot(Pf_wac       ,
            x = 0.01 , y = .88, width = .19, height = .11)+
  draw_plot(P_wac         ,
            x = .0, y = .66, width = .2, height = .22)+
  
  draw_plot(Pf_block       ,
            x = .25 ,y = .88, width = .15, height = .11)+
  draw_plot(P_block         ,
            x = .2, y = .66, width = .25, height = .22)+
  draw_plot(Pf_decay       ,
            x = .45 , y = .88, width = .15, height = .11)+
  draw_plot(P_decay         ,
            x = .4, y = .66, width = .25, height = .22)+
  
   
  
  draw_plot(P1         + theme(legend.position = "none"),
            x = 0.0 , y = .33, width = .2, height = .33)+
  draw_plot(P2         + theme(legend.position = "none"),
            x = .2, y = .33, width = .2, height = .33)+
  draw_plot(P3         + theme(legend.position = "none"),
            x = .4, y = .33, width = .2, height = .33)+
  draw_plot(P4         + theme(legend.position = "none"),
            x = 0 , y = .0, width = .2, height = .33)+
  draw_plot(P5         + theme(legend.position = "none"),
            x = .2, y = .0, width = .2, height = .33)+
  draw_plot(P6         + theme(legend.position = "none"),
            x = .4, y = .0, width = .2, height = .33)+
  draw_plot(P1_p       + theme(legend.position = "none"),
            x = .6, y = .6, width = .2, height = .35)+
  draw_plot(P1_t1      + theme(legend.position = "none"),
            x = .8, y = .6, width = .2, height = .35)+
  draw_plot(p_run_time + theme(legend.position = "none"),
            x = .6, y = .1, width = .4, height = .4)+
  
  draw_label("A", x = 0.01, y = 0.98, vjust = 1  ) +
  
  draw_label("B", x = 0.01, y = 0.62, vjust = 1  ) +
  
  draw_label("D", x =0.62, y = 0.98,  vjust = 1) +
  
  draw_label("C", x = 0.01, y = .33,   vjust = 1)+
  draw_label("E", x = 0.62, y = .5,   vjust = 1)
legend <- get_legend(
  P1 +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)
P_out <- plot_grid(grid_plot, legend , ncol=1,rel_heights = c(0.9, .1))


ggsave(P_out , file="plot/Fig3.png",
       width = 29.7,
       height = 21,
       units = "cm"
)
