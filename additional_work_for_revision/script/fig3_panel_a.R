library(dplR)
library(ggplot2)
library(gridExtra)
library(grid)

library(cowplot)
library(dplyr)
library(susieR)
library(ggpubr)
#D41159  
#1A85FF
#40B0A6

colors <- c("gold","#D41159","dodgerblue" )
### For PVE=10% ----

simple_roc <- function(labs, scores){
  labs <- labs[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labs)/sum(labs), FPR=cumsum(!labs)/sum(!labs), labs)
}

#### preparing the data  ----
path <- getwd()
load(paste( path,"/simulation/Simulation_results/comparison_susie_fusie_128_sd1.RData", 
            sep=""))


run_time_fsusie <- do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_time[3]))

run_time_sp_fsusie <- do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_sp_time[3]))

run_time_susie <- do.call( c, lapply( 1:length(res), function( i) res[[i]]$susie_time[3]))




idx <- 1:length(res)
#idx <- 20:50
true_lab <- do.call( c,
                     lapply(idx,
                            
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


score_susie <-  do.call( c, lapply(idx,
                                   function( i) res[[i]]$susie_rss_pip))
score_fsusie <-  do.call( c, lapply( idx,
                                     function( i) res[[i]]$susiF_pip))
score_sp_fsusie <-  do.call( c, lapply(idx,
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
 


#### ROC Gaussian simulation   ------

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
  geom_line(linewidth=1.2)+
  xlim( c(0,0.05))+
  theme(legend.position = "none")+
  theme_linedraw()+
  ggtitle("Wavelet")+
  scale_color_manual(values = colors)
P1

#### ROCFractional   ------#long_range_128_sd3_pip.RData")
load("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/additional_work_for_revision/simulations/pip_res/long_range_128_sd3_pip.RData")
true_lab <- do.call( c,
                     lapply(1: length(res),
                            
                            function( i) {
                              
                              a <-  rep( 0,   length(res[[i]]$pip_fsusie_ps))
                              a[res[[i]]$true_pos] <- 1
                              return(a)
                            }
                            
                     )
)  


simple_roc <- function(labs, scores){
  labs <- labs[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labs)/sum(labs), FPR=cumsum(!labs)/sum(!labs), labs)
}
#$pip_susie
#$pip_fsusie_ps
#$pip_fsusie_nps

score_fsusie <-  do.call( c, lapply( 1: length(res),
                                     function( i) res[[i]]$pip_fsusie_nps))
score_sp_fsusie <-  do.call( c, lapply( 1: length(res),
                                        function( i) res[[i]]$pip_fsusie_ps))

score_susie <-  do.call( c, lapply( 1: length(res),
                                    function( i) res[[i]]$pip_susie))



pip_scores <- score_fsusie 
# 
causal_flags <-true_lab
 

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

P2 <- ggplot(df_roc, aes (x=FDR, y=Power,col=method))+
  geom_line(linewidth=1.2)+
  xlim( c(0,0.05))+
  theme(legend.position = "none")+
  theme_linedraw()+
  ggtitle("Wavelet fractional noise")+
  scale_color_manual(values = colors)
P2


#### ROCGP   ------
load("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/additional_work_for_revision/simulations/pip_res/short_range_128_sd3_pip.RData")
true_lab <- do.call( c,
                     lapply(1: length(res),
                            
                            function( i) {
                              
                              a <-  rep( 0,   length(res[[i]]$pip_fsusie_ps))
                              a[res[[i]]$true_pos] <- 1
                              return(a)
                            }
                            
                     )
)  


simple_roc <- function(labs, scores){
  labs <- labs[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labs)/sum(labs), FPR=cumsum(!labs)/sum(!labs), labs)
}
#$pip_susie
#$pip_fsusie_ps
#$pip_fsusie_nps

score_fsusie <-  do.call( c, lapply( 1: length(res),
                                     function( i) res[[i]]$pip_fsusie_nps))
score_sp_fsusie <-  do.call( c, lapply( 1: length(res),
                                        function( i) res[[i]]$pip_fsusie_ps))

score_susie <-  do.call( c, lapply( 1: length(res),
                                    function( i) res[[i]]$pip_susie))



pip_scores <- score_fsusie 
# 
causal_flags <-true_lab


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

P3 <- ggplot(df_roc, aes (x=FDR, y=Power,col=method))+
  geom_line(linewidth=1.2)+
  xlim( c(0,0.05))+
  theme(legend.position = "none")+
  theme_linedraw()+
  ggtitle("Wavelet GP noise")+
  scale_color_manual(values = colors)
P3

#### ROC block simulation  ----
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
score_sp_fsusie <-  do.call( c, lapply( 1: length(res),
                                        function( i) res[[i]]$susiF_sp_pip))


 
 
  
 
score_susie <-  do.call( c, lapply( 1: length(res),
                                    function( i) res[[i]]$susie_rss_pip))




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
P4 <-ggplot(df_roc, aes (x=FDR, y=Power,col=method))+
  geom_line(linewidth=1.2)+
  xlim( c(0,0.05))+
  theme(legend.position = "none")+
  theme_linedraw()+
  ggtitle("WGBS block")+
  scale_color_manual(values = colors)


P4
#### ROC  distance decay simulation  -----

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
data(N3finemapping)
X <- N3finemapping$X
True_cor <- cor(X)







score_fsusie <-  do.call( c, lapply( 1: length(res),
                                     function( i) res[[i]]$susiF_pip))
score_sp_fsusie <-  do.call( c, lapply( 1: length(res),
                                        function( i) res[[i]]$susiF_sp_pip))


 

score_susie <-  do.call( c, lapply( 1: length(res),
                                    function( i) res[[i]]$susie_rss_pip))



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
P5 <-ggplot(df_roc, aes (x=FDR, y=Power,col=method))+
  geom_line(linewidth=1.2)+
  xlim( c(0,0.05))+
  theme_linedraw()+
  ggtitle("WGBS distance decay")+
  scale_color_manual(values = colors)
P5

 

ggarrange(P1,P2,P3,P4,P5, ncol=5, common.legend = TRUE, legend="bottom")

save_path ="C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/additional_work_for_revision/plot/"
ggsave(ggarrange(P1,P2,P3,P4,P5, ncol=5, common.legend = TRUE, legend="bottom"), 
       file = paste0(save_path, "power.pdf"),
       width = 30.5, height = 10.5, units = "cm")
