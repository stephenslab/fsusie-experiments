load("D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/simulation/Simulation_script/script/ROC_sim_CpG_discovery/h2_05.RData")
library( ggplot2)
library(cowplot)
library(gridExtra)
res= do.call( rbind, res)



simple_roc <- function(labs, scores){
  labs <- labs[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labs)/sum(labs), FPR=cumsum(!labs)/sum(!labs), labs)
}

#usage 
roc0 <- simple_roc(res$CpG,   res$pv)
plot( roc$TPR,roc$FPR, type="l", lwd=2)
labs  <- res$CpG
score <-   res$hmm_lfsr
roc1 <- simple_roc(res$CpG, res$hmm_lfsr)

lines( roc$TPR,roc$FPR, col="blue", lwd=2)
labs  <- res$CpG
score <-   res$pv_ti
roc2 <- simple_roc(res$CpG,res$pv_ti)
df_plot = data.frame (TPR= c(roc0$TPR,
                              roc1$TPR,
                              roc2$TPR),
                      FPR=  c(roc0$FPR,
                              roc1$FPR,
                              roc2$FPR),
                      col=factor(c( rep( "pv", length(roc0$FPR )),
                                    rep( "hmm", length(roc1$FPR )),
                                    rep( "TI", length(roc1$FPR )))
                                )
)

ggplot(df_plot, aes(x=TPR, y=FPR, col=col))+geom_line()+
  theme_cowplot()+
  ggtitle()



lines( roc$TPR,roc$FPR, col="lightgreen", lwd=2)

load("D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/simulation/Simulation_script/script/ROC_sim_CpG_discovery/h2_10.RData")

res= do.call( rbind, res)



simple_roc <- function(labs, scores){
  labs <- labs[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labs)/sum(labs), FPR=cumsum(!labs)/sum(!labs), labs)
}

#usage
labs  <- res$CpG
score <-   res$pv
roc <- simple_roc(labs, score)
plot( roc$TPR,roc$FPR, type="l", lwd=2)
labs  <- res$CpG
score <-   res$hmm_lfsr
roc <- simple_roc(labs, score)

lines( roc$TPR,roc$FPR, col="blue", lwd=2)
labs  <- res$CpG
score <-   res$pv_ti
roc <- simple_roc(labs, score)

lines( roc$TPR,roc$FPR, col="lightgreen", lwd=2)





load("D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/simulation/Simulation_script/script/ROC_sim_CpG_discovery/h2_20.RData")

res= do.call( rbind, res)



simple_roc <- function(labs, scores){
  labs <- labs[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labs)/sum(labs), FPR=cumsum(!labs)/sum(!labs), labs)
}

#usage
labs  <- res$CpG
score <-   res$pv
roc <- simple_roc(labs, score)
plot( roc$TPR,roc$FPR, type="l", lwd=2)
labs  <- res$CpG
score <-   res$hmm_lfsr
roc <- simple_roc(labs, score)

lines( roc$TPR,roc$FPR, col="blue", lwd=2)
labs  <- res$CpG
score <-   res$pv_ti
roc <- simple_roc(labs, score)

lines( roc$TPR,roc$FPR, col="lightgreen", lwd=2)

