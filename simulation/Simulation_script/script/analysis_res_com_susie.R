library(dplR)
library(ggplot2)
library(gridExtra)
library(grid)


library(dplyr)
library(susieR)
library(ggpubr)
### For PVE=10% ----
load("~/fsusi_simu/sim3/comparison_susie_fusie_128_sd1.RData")


run_time_fsusie <- do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_time[3]))

run_time_sp_fsusie <- do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_sp_time[3]))

run_time_susie <- do.call( c, lapply( 1:length(res), function( i) res[[i]]$susie_time[3]))




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
score_sp_fsusie <-  do.call( c, lapply( 1: length(res),
                                        function( i) res[[i]]$susiF_sp_pip))

cs_size_fsusie1 <- mean(  do.call( c, lapply( 1: length(res),
                                              function( i)   lengths(res[[i]]$susiF_cs)) ))

cs_size_sp_fsusie1 <- mean(  do.call( c, lapply( 1: length(res),
                                                 function( i)   lengths(res[[i]]$susiF_sp_cs)) ))


cs_size_susie1 <- mean(  do.call( c, lapply( 1: length(res),
                                             function (i){
                                               
                                               if  (length(res[[i]]$susie_cs$cs ) ==0 ){
                                                 return(NA)
                                               }else{
                                                 out <-  do.call (c, lapply(1:length(res[[i]]$susie_cs$cs) ,
                                                                            function(l ) length(res[[i]]$susie_cs$cs[[l]]) ))
                                               }
                                               return( out)
                                             }
                                             
)
), na.rm=TRUE
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


tt <- do.call( c, lapply( 1: length(res),
                          function( i)  lengths(res[[i]]$susiF_sp_cs)))
t0 <-  do.call( c, lapply( 1: length(res),
                           function( i) res[[i]]$fsusie_sp_purity))
if (  length(which(t0<0.60))>0){
  
  purity_fsusie_sp_1  <-mean( t0[-which(t0<0.60)])
  
}else{
  
  purity_fsusie_sp_1  <-mean( t0 )
  
}

if (  length(which(tt>50))>0){
  
  cs_size_fsusie_sp_1  <-mean( tt[-which(tt>50)])
  
}else{
  
  cs_size_fsusie_sp_1   <-mean( tt )
  
}




n_effect <-  do.call(c, lapply( 1: length(res),
                                function( i) length(res[[i]]$ true_pos)))

score_susie <-  do.call( c, lapply( 1: length(res),
                                    function( i) res[[i]]$susie_rss_pip))

simple_roc <- function(labs, scores){
  labs <- labs[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labs)/sum(labs), FPR=cumsum(!labs)/sum(!labs), labs)
}

## corriger sous ça ------

roc_fsusie <- simple_roc(true_lab, score_fsusie)
roc_sp_fsusie <- simple_roc(true_lab, score_sp_fsusie  )
roc_susie <- simple_roc(true_lab, score_susie)



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
  ggtitle("Gaussian functional effect ")
P1






### For PVE=10% ----
load("~/fsusi_simu/sim3/comparison_susie_fusie_block_sd1.RData")

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




run_time_fsusie <- c(run_time_fsusie ,do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_time[3])) )

run_time_sp_fsusie <-c(run_time_sp_fsusie,  do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_sp_time[3])))

run_time_susie <- c(run_time_susie ,do.call( c, lapply( 1:length(res), function( i) res[[i]]$susie_time[3])))





score_fsusie <-  do.call( c, lapply( 1: length(res),
                                     function( i) res[[i]]$susiF_pip))
score_sp_fsusie <-  do.call( c, lapply( 1: length(res),
                                        function( i) res[[i]]$susiF_sp_pip))

cs_size_fsusie1 <- mean(  do.call( c, lapply( 1: length(res),
                                              function( i)   lengths(res[[i]]$susiF_cs)) ))

cs_size_sp_fsusie1 <- mean(  do.call( c, lapply( 1: length(res),
                                                 function( i)   lengths(res[[i]]$susiF_sp_cs)) ))


cs_size_susie2 <- mean(  do.call( c, lapply( 1: length(res),
                                             function (i){
                                               
                                               if  (length(res[[i]]$susie_cs$cs ) ==0 ){
                                                 return(NA)
                                               }else{
                                                 out <-  do.call (c, lapply(1:length(res[[i]]$susie_cs$cs) ,
                                                                            function(l ) length(res[[i]]$susie_cs$cs[[l]]) ))
                                               }
                                               return( out)
                                             }
                                             
)
), na.rm=TRUE
)



purity_susie2 <- mean(  do.call( c, lapply( 1: length(res),
                                            function( i)  res[[i]]$susie_cs$purity [,1]  )), na.rm = TRUE)

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


tt <- do.call( c, lapply( 1: length(res),
                          function( i)  lengths(res[[i]]$susiF_sp_cs)))
t0 <-  do.call( c, lapply( 1: length(res),
                           function( i) res[[i]]$fsusie_sp_purity))
if (  length(which(t0<0.60))>0){
  
  purity_fsusie_sp_2  <-mean( t0[-which(t0<0.60)])
  
}else{
  
  purity_fsusie_sp_2  <-mean( t0 )
  
}

if (  length(which(tt>50))>0){
  
  cs_size_fsusie_sp_2  <-mean( tt[-which(tt>50)])
  
}else{
  
  cs_size_fsusie_sp_2   <-mean( tt )
  
}




n_effect <-  do.call(c, lapply( 1: length(res),
                                function( i) length(res[[i]]$ true_pos)))

score_susie <-  do.call( c, lapply( 1: length(res),
                                    function( i) res[[i]]$susie_rss_pip))

simple_roc <- function(labs, scores){
  labs <- labs[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labs)/sum(labs), FPR=cumsum(!labs)/sum(!labs), labs)
}


roc_fsusie <- simple_roc(true_lab, score_fsusie)
roc_sp_fsusie <- simple_roc(true_lab, score_sp_fsusie  )
roc_susie <- simple_roc(true_lab, score_susie)



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
  ggtitle("WGBS block effect ")


P2
### For PVE=10% ----
load("~/fsusi_simu/sim3/comparison_susie_fusie_distdecay_sd1.RData")

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


run_time_fsusie <- c(run_time_fsusie ,do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_time[3])) )

run_time_sp_fsusie <-c(run_time_sp_fsusie,  do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_sp_time[3])))

run_time_susie <- c(run_time_susie ,do.call( c, lapply( 1:length(res), function( i) res[[i]]$susie_time[3])))




score_fsusie <-  do.call( c, lapply( 1: length(res),
                                     function( i) res[[i]]$susiF_pip))
score_sp_fsusie <-  do.call( c, lapply( 1: length(res),
                                        function( i) res[[i]]$susiF_sp_pip))

cs_size_fsusie1 <- mean(  do.call( c, lapply( 1: length(res),
                                              function( i)   lengths(res[[i]]$susiF_cs)) ))

cs_size_sp_fsusie1 <- mean(  do.call( c, lapply( 1: length(res),
                                                 function( i)   lengths(res[[i]]$susiF_sp_cs)) ))


cs_size_susie3 <- mean(  do.call( c, lapply( 1: length(res),
                                             function (i){
                                               
                                               if  (length(res[[i]]$susie_cs$cs ) ==0 ){
                                                 return(NA)
                                               }else{
                                                 out <-  do.call (c, lapply(1:length(res[[i]]$susie_cs$cs) ,
                                                                            function(l ) length(res[[i]]$susie_cs$cs[[l]]) ))
                                               }
                                               return( out)
                                             }
                                             
)
), na.rm=TRUE
)



purity_susie3 <- mean(  do.call( c, lapply( 1: length(res),
                                            function( i)  res[[i]]$susie_cs$purity [,1]  )), na.rm = TRUE)

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


tt <- do.call( c, lapply( 1: length(res),
                          function( i)  lengths(res[[i]]$susiF_sp_cs)))
t0 <-  do.call( c, lapply( 1: length(res),
                           function( i) res[[i]]$fsusie_sp_purity))
if (  length(which(t0<0.60))>0){
  
  purity_fsusie_sp_3  <-mean( t0[-which(t0<0.60)])
  
}else{
  
  purity_fsusie_sp_3  <-mean( t0 )
  
}

if (  length(which(tt>50))>0){
  
  cs_size_fsusie_sp_3  <-mean( tt[-which(tt>50)])
  
}else{
  
  cs_size_fsusie_sp_3   <-mean( tt )
  
}




n_effect <-  do.call(c, lapply( 1: length(res),
                                function( i) length(res[[i]]$ true_pos)))

score_susie <-  do.call( c, lapply( 1: length(res),
                                    function( i) res[[i]]$susie_rss_pip))

simple_roc <- function(labs, scores){
  labs <- labs[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labs)/sum(labs), FPR=cumsum(!labs)/sum(!labs), labs)
}

## corriger sous ça ------

roc_fsusie <- simple_roc(true_lab, score_fsusie)
roc_sp_fsusie <- simple_roc(true_lab, score_sp_fsusie  )
roc_susie <- simple_roc(true_lab, score_susie)



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
  ggtitle("WGBS distance decay effect ")
P3

df_cs_purity <- data.frame( scenario= factor(
  rep(
    c("Gaussian","WGBS block", "WGBS decay"),
    each=3)),
  cs_size=c(cs_size_fsusie1, cs_size_susie1, cs_size_sp_fsusie1,
            cs_size_fsusie2, cs_size_susie2, cs_size_fsusie2,
            cs_size_fsusie3, cs_size_susie3, cs_size_fsusie3),
  purity=c(purity_fsusie1, purity_susie1,purity_fsusie_sp_1,
           purity_fsusie2, purity_susie2,purity_fsusie_sp_2,
           purity_fsusie3, purity_susie3,purity_fsusie_sp_3),
  method= factor( rep(c("fSuSiE SPS","SuSiE","fSuSiE IS"),3))
)
P4 <- ggplot( df_cs_purity, aes(x= scenario, y=cs_size, col=method))+
  geom_point(size=2)+ scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_linedraw() +
  xlab("")+
  ylab("CS size")
P5  <-ggplot( df_cs_purity, aes(x= scenario, y=purity, col=method))+
  geom_point(size=2)+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_linedraw() +
  xlab("")+
  ylab("Purity")





## run time df -----
df_run_time = data.frame( run_time=  c(run_time_fsusie,
                                       run_time_sp_fsusie,
                                       128*run_time_susie),
                          method=  factor ( c(rep( "fSuSiE SPS", length( run_time_fsusie)),
                                              rep( "fSuSiE IS", length(run_time_sp_fsusie)),
                                              rep( "SuSiE", length( run_time_susie))
                          )
                          )
)

P_run_time <- ggplot( df_run_time, aes( x= run_time, y=method, col=method))+
  geom_boxplot()+ theme_linedraw() +ylab("method")+xlab("Run time")



load("~/fsusi_simu/sim3/overlap_check_gaussian.RData")

mean_overlapp_susif_gaus <- mean(do.call( c, lapply(1: length(res), function(i) res[[i]]$is_overlap_susif)))
mean_overlapp_susie_gaus <-mean( do.call( c, lapply(1: length(res), function(i) res[[i]]$is_overlap_susie)))

mean_overlapp_susif_sp_gaus <- mean(do.call( c, lapply(1: length(res), function(i) res[[i]]$is_overlap_susif_sp)))

load("~/fsusi_simu/sim3/overlap_check_block_sd1.RData")

mean_overlapp_susif_block <- mean(do.call( c, lapply(1: length(res), function(i) res[[i]]$is_overlap_susif)))

mean_overlapp_susie_block  <-mean( do.call( c, lapply(1: length(res), function(i) res[[i]]$is_overlap_susie)))
mean_overlapp_susif_sp_block <- mean(do.call( c, lapply(1: length(res), function(i) res[[i]]$is_overlap_susif_sp)))

load("~/fsusi_simu/sim3/overlap_check_distdecay_sd1.RData")

mean_overlapp_susif_decay <- mean(do.call( c, lapply(1: length(res), function(i) res[[i]]$is_overlap_susif)))
mean_overlapp_susie_decay  <-mean( do.call( c, lapply(1: length(res), function(i) res[[i]]$is_overlap_susie)))
mean_overlapp_susif_sp_decay <- mean(do.call( c, lapply(1: length(res), function(i) res[[i]]$is_overlap_susif_sp)))

df_overlapp <- data.frame( scenario= factor(
  rep(
    c("Gaussian","WGBS block", "WGBS decay"),
    each=3)),
  overlapp=c(mean_overlapp_susif_gaus ,mean_overlapp_susie_gaus,mean_overlapp_susif_sp_gaus ,
             mean_overlapp_susif_block, mean_overlapp_susie_block ,mean_overlapp_susif_sp_block ,
             mean_overlapp_susif_decay,   mean_overlapp_susie_decay,mean_overlapp_susif_sp_decay ),
  
  method= factor( rep(c("fSuSiE SPS","SuSiE","fSuSiE IS"),3))
)
P6  <-ggplot(df_overlapp, aes(x= scenario, y=overlapp, col=method))+
  geom_point(size=2)+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_linedraw() +
  xlab("")+
  ylab("Probability of overlapp")+scale_y_continuous(labels = scales::percent)




#####


library(dplR)
### For PVE=10% ----
load("/home/wdenault/fsusi_simu/sim3/check_L_accuracy_128_sd1.RData")




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




P_out <-ggarrange(P1,P2,P3,P4,P5,P6 , P1_p, P1_t1, ncol=4, nrow=2, common.legend = TRUE, legend="bottom")
P_out




##### run time ------
load("~/fsusi_simu/sim3/run_time_comp_p100.RData")


run_time_fsusie100 <- do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_time[3]))

run_time_sp_fsusie100  <- do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_sp_time[3]))


load("~/fsusi_simu/sim3/run_time_comp_p500.RData")


run_time_fsusie500 <- do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_time[3]))

run_time_sp_fsusie500  <- do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_sp_time[3]))


load("~/fsusi_simu/sim3/run_time_comp_p1000.RData")


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
p_run_time <- ggplot(df_run_time, aes(runtime ,y= prior, col=prior))+
  geom_boxplot()+
  facet_wrap( N~.,scale="free")+
  theme_linedraw() +
  ylab("Prior")+
  xlab("Run time")+
  theme(strip.background =element_rect(fill="white"),
        strip.text = element_text(colour = 'black'))






grid_plot <- ggdraw()+
  draw_plot(P1         + theme(legend.position = "none"),
            x = 0 , y = .5, width = .2, height = .5)+
  draw_plot(P2         + theme(legend.position = "none"),
            x = .2, y = .5, width = .2, height = .5)+
  draw_plot(P3         + theme(legend.position = "none"),
            x = .4, y = .5, width = .2, height = .5)+
  draw_plot(P4         + theme(legend.position = "none"),
            x = 0 , y = .0, width = .2, height = .5)+
  draw_plot(P5         + theme(legend.position = "none"),
            x = .2, y = .0, width = .2, height = .5)+
  draw_plot(P6         + theme(legend.position = "none"),
            x = .4, y = .0, width = .2, height = .5)+
  draw_plot(P1_p       + theme(legend.position = "none"),
            x = .6, y = .5, width = .2, height = .5)+
  draw_plot(P1_t1      + theme(legend.position = "none"),
            x = .8, y = .5, width = .2, height = .5)+
  draw_plot(p_run_time + theme(legend.position = "none"),
            x = .6, y = .0, width = .4, height = .5)
legend <- get_legend(
  P1 +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)
plot_grid(grid_plot, legend , ncol=1,rel_heights = c(0.9, .1))
