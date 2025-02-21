library(dplR)
library(ggplot2)
library(gridExtra)
library(grid)


library(dplyr)
library(susieR)
library(ggpubr)
#D41159  
#1A85FF
#40B0A6

colors <- c("#D41159","#1A85FF","#40B0A6" )
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
  geom_line(size=1.2)+
  xlim( c(0,0.05))+
  theme(legend.position = "none")+
  theme_linedraw()+
  ggtitle("Gaussian functional")+
  scale_color_manual(values = colors)
P1



### For PVE=10% ----

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
  ggtitle("WGBS block")+
  scale_color_manual(values = colors)


P2
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
  ggtitle("WGBS distance decay")+
  scale_color_manual(values = colors)
P3



#### purity and CS size plots-----

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



#### run time plot -----
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



#### Overlapp probability plot ----

path <- getwd()
load(paste( path,"/simulation/Simulation_results/overlap_check_distdecay_sd1.RData", 
            sep=""))


mean_overlapp_susif_gaus <- mean(do.call( c, lapply(1: length(res), function(i) res[[i]]$is_overlap_susif)))
mean_overlapp_susie_gaus <-mean( do.call( c, lapply(1: length(res), function(i) res[[i]]$is_overlap_susie)))

mean_overlapp_susif_sp_gaus <- mean(do.call( c, lapply(1: length(res), function(i) res[[i]]$is_overlap_susif_sp)))

path <- getwd()
load(paste( path,"/simulation/Simulation_results/overlap_check_block_sd1.RData", 
            sep=""))


mean_overlapp_susif_block <- mean(do.call( c, lapply(1: length(res), function(i) res[[i]]$is_overlap_susif)))

mean_overlapp_susie_block  <-mean( do.call( c, lapply(1: length(res), function(i) res[[i]]$is_overlap_susie)))
mean_overlapp_susif_sp_block <- mean(do.call( c, lapply(1: length(res), function(i) res[[i]]$is_overlap_susif_sp)))


path <- getwd()
load(paste( path,"/simulation/Simulation_results/overlap_check_distdecay_sd1.RData", 
            sep=""))

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
  geom_point(size=3)+
  scale_x_discrete(guide = guide_axis(angle = 30)) +
  theme_linedraw() +
  xlab("")+
  
  ylab("Probability of overlapp")+
  scale_y_continuous(labels = scales::percent,
                     limits =c(0.65,1.0))+
  scale_color_manual(values = colors)


##### ROC plot ----

grid_plot_roc <- ggdraw()+

draw_plot(P1         + theme(legend.position = "none", 
                             title = element_text( size=10 )),
          x = 0.0 , y = .0, width = .33, height = 1)+
  draw_plot(P2         + theme(legend.position = "none", 
                               title = element_text( size=10 )),
            x = .33, y = .0, width = .33, height =1)+
  draw_plot(P3         + theme(legend.position = "none", 
                               title = element_text( size=10 )),
            x = 0.66 , y = .0, width = .33, height = 1)

 
save_path=  paste0(getwd(),
                   "/plot/fig3_separate_panel/"
)
ggsave(grid_plot_roc , file=paste0(save_path,"ROC.pdf"),
       width = 29.7,
       height = 21,
       units = "cm"
)


grid_plot_summary <- ggdraw()+
  
  draw_plot(P4         + theme(legend.position = "none", 
                               title = element_text( size=10 )),
            x = 0.0 , y = .0, width = .33, height = 1)+
  draw_plot(P5         + theme(legend.position = "none", 
                               title = element_text( size=10 )),
            x = .33, y = .0, width = .33, height =1)+
  draw_plot(P6         + theme(legend.position = "none", 
                               title = element_text( size=10 )),
            x = 0.66 , y = .0, width = .33, height = 1)



save_path=  paste0(getwd(),
                   "/plot/fig3_separate_panel/"
)
ggsave(grid_plot_summary , file=paste0(save_path,"purity_cs_size_overlapp.pdf"),
       width = 29.7,
       height = 21,
       units = "cm"
)
