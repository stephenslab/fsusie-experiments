library(dplR)
### For PVE=10% ----
load("~/fsusi_simu/sim3/comparison_susie_fusie_128_sd1.RData")
true_lab <- do.call( c,
                     lapply(1: length(res),

                     function( i) {

                       a <-  rep( 0,   length(res[[1]]$susiF_pip))
                       a[res[[i]]$true_pos] <- 1
                       return(a)
                     }

                     )
)



score_fsusie <-  do.call( c, lapply( 1: length(res),
                                     function( i) res[[i]]$susiF_pip))

score_susie <-  do.call( c, lapply( 1: length(res),
                                     function( i) res[[i]]$susie_rss_pip))

 simple_roc <- function(labs, scores){
   labs <- labs[order(scores, decreasing=TRUE)]
   data.frame(TPR=cumsum(labs)/sum(labs), FPR=cumsum(!labs)/sum(!labs), labs)
 }



 roc_fsusie <- simple_roc(true_lab, score_fsusie)
 roc_susie <- simple_roc(true_lab, score_susie)
df_roc <- data.frame (FDR  =c( roc_fsusie$FPR, roc_susie$FPR),
                      Power = c( roc_fsusie$TPR, roc_susie$TPR),
                      method= factor ( c(rep("fSuSIE", length(roc_fsusie$FPR)),
                                         rep("SuSIE", length(roc_susie$FPR))
                                         )
                                       )
                      )
library(ggplot2)
ggplot(df_roc, aes (x=FDR, y=Power,col=method))+
  geom_line()+
xlim( c(0,0.2))



library(dplR)
### For PVE=10% ----
load("~/fsusi_simu/sim3/comparison_susie_fusie_distdecay_sd1.RData")
true_lab <- do.call( c,
                     lapply(1: length(res),

                            function( i) {

                              a <-  rep( 0,   length(res[[1]]$susiF_pip))
                              a[res[[i]][[5]]] <- 1
                              return(a)
                            }

                     )
)



score_fsusie <-  do.call( c, lapply( 1: length(res),
                                     function( i) res[[i]]$susiF_pip))

score_susie <-  do.call( c, lapply( 1: length(res),
                                    function( i) res[[i]]$susie_rss_pip))


scores <-score_fsusie
labs  <- true_lab
simple_roc <- function(labs, scores){
  labs <- labs[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labs)/sum(labs), FPR=cumsum(!labs)/sum(!labs), labs)
}



roc_fsusie <- simple_roc(true_lab, score_fsusie)
roc_susie <- simple_roc(true_lab, score_susie)
df_roc <- data.frame (FDR  =c( roc_fsusie$FPR, roc_susie$FPR),
                      Power = c( roc_fsusie$TPR, roc_susie$TPR),
                      method= factor ( c(rep("fSuSIE", length(roc_fsusie$FPR)),
                                         rep("SuSIE", length(roc_susie$FPR))
                      )
                      )
)
library(ggplot2)
ggplot(df_roc, aes (x=FDR, y=Power,col=method))+
  geom_line()+
  xlim( c(0,0.2))




library(dplR)
### For PVE=10% ----
load("~/fsusi_simu/sim3/comparison_susie_fusie_distdecay_sd1.RData")
true_lab <- do.call( c,
                     lapply(1: length(res),

                            function( i) {

                              a <-  rep( 0,   length(res[[1]]$susiF_pip))
                              a[res[[i]][[5]]] <- 1
                              return(a)
                            }

                     )
)



score_fsusie <-  do.call( c, lapply( 1: length(res),
                                     function( i) res[[i]]$susiF_pip))

score_susie <-  do.call( c, lapply( 1: length(res),
                                    function( i) res[[i]]$susie_rss_pip))


scores <-score_fsusie
labs  <- true_lab
simple_roc <- function(labs, scores){
  labs <- labs[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labs)/sum(labs), FPR=cumsum(!labs)/sum(!labs), labs)
}



roc_fsusie <- simple_roc(true_lab, score_fsusie)
roc_susie <- simple_roc(true_lab, score_susie)
df_roc <- data.frame (FDR  =c( roc_fsusie$FPR, roc_susie$FPR),
                      Power = c( roc_fsusie$TPR, roc_susie$TPR),
                      method= factor ( c(rep("fSuSIE", length(roc_fsusie$FPR)),
                                         rep("SuSIE", length(roc_susie$FPR))
                      )
                      )
)
library(ggplot2)
ggplot(df_roc, aes (x=FDR, y=Power,col=method))+
  geom_line()+
  xlim( c(0,0.2))

