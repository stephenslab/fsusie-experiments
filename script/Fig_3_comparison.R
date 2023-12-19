library(dplR)
library(ggplot2)
library(gridExtra)
library(grid)


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



roc_fsusie <- simple_roc(true_lab, score_fsusie)
roc_susie <- simple_roc(true_lab, score_susie)
df_roc <- data.frame ( Power =c( roc_fsusie$TPR, roc_susie$TPR),
                       FDR = c( roc_fsusie$FPR, roc_susie$FPR),
                       method= factor ( c(rep("fSuSIE", length(roc_fsusie$FPR)),
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



roc_fsusie <- simple_roc(true_lab, score_fsusie)
roc_susie <- simple_roc(true_lab, score_susie)
df_roc <- data.frame (Power  =c( roc_fsusie$TPR, roc_susie$TPR),
                      FDR = c( roc_fsusie$FPR, roc_susie$FPR),
                      method= factor ( c(rep("fSuSIE", length(roc_fsusie$FPR)),
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





scores <-score_fsusie
labs  <- true_lab
simple_roc <- function(labs, scores){
  labs <- labs[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labs)/sum(labs), FPR=cumsum(!labs)/sum(!labs), labs)
}



roc_fsusie <- simple_roc(true_lab, score_fsusie)
roc_susie <- simple_roc(true_lab, score_susie)
df_roc <- data.frame (  Power =c( roc_fsusie$TPR, roc_susie$TPR),
                        FDR = c( roc_fsusie$FPR, roc_susie$FPR),
                        method= factor ( c(rep("fSuSIE", length(roc_fsusie$FPR)),
                                           rep("SuSIE", length(roc_susie$FPR))
                        )
                        )
)
P3 <-ggplot(df_roc, aes (x=FDR, y=Power,col=method))+
  geom_line(size=1.2)+
  xlim( c(0,0.05))+
  theme_linedraw()+
  ggtitle("WGBS distance decay effect ")


df_cs_purity <- data.frame( scenario= factor(
  rep(
    c("Gaussian","WGBS block", "WGBS decay"),
    each=2)),
  cs_size=c(cs_size_fsusie1, cs_size_susie1,
            cs_size_fsusie2, cs_size_susie2,
            cs_size_fsusie3, cs_size_susie3),
  purity=c(purity_fsusie1, purity_susie1,
           purity_fsusie2, purity_susie2,
           purity_fsusie3, purity_susie3),
  method= factor( rep(c("fSuSiE","SuSiE"),3))
)
P4 <- ggplot( df_cs_purity, aes(x= scenario, y=cs_size, col=method))+
  geom_point(size=3)+ scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_linedraw() +
  xlab("")+
  ylab("CS size")
P5  <-ggplot( df_cs_purity, aes(x= scenario, y=purity, col=method))+
  geom_point(size=3)+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_linedraw() +
  xlab("")+
  ylab("Purity")



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
    each=2)),
  overlapp=c(mean_overlapp_susif_gaus ,mean_overlapp_susie_gaus,
             mean_overlapp_susif_block, mean_overlapp_susie_block ,
             mean_overlapp_susif_decay,   mean_overlapp_susie_decay),
  
  method= factor( rep(c("fSuSiE","SuSiE"),3))
)
P6  <-ggplot(df_overlapp, aes(x= scenario, y=overlapp, col=method))+
  geom_point(size=3)+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_linedraw() +
  xlab("")+
 
  ylab("Probability of overlapp")+
  scale_y_continuous(labels = scales::percent,
                     limits =c(0.65,1.0))




P_out <- ggarrange(
                   P1,P2,P3,P4,P5,P6, 
                   ncol=3, 
                   nrow=2, 
                   common.legend = TRUE, 
                   legend="bottom"
                   )
P_out
ggsave(P_out , file="plot/Fig3.png",
       width = 29.7,
       height = 21,
       units = "cm"
)
