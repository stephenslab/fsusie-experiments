load("D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/simulation/Simulation_script/script/ROC_sim_CpG_discovery/h2_01_n100_CpG_10_decay.RData")
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
roc1 <- simple_roc(res$CpG, res$hmm_lfsr)
roc2 <- simple_roc(res$CpG,res$pv_ti)
roc3 <- simple_roc(res$CpG,res$affected_TI)
df_plot = data.frame (TPR= c(roc0$TPR,
                             roc1$TPR,
                             roc2$TPR ),
                      FPR=  c(roc0$FPR,
                              roc1$FPR,
                              roc2$FPR ),
                      col=factor(c( rep( "pv", length(roc0$FPR )),
                                    rep( "hmm", length(roc1$FPR )),
                                    rep( "TI_pv", length(roc1$FPR )) )
                      )
)

P11= ggplot(df_plot, aes(x=TPR, y=FPR, col=col))+geom_line()+
  theme_cowplot()+xlab("FDR")+ylab("Power")+
  ggtitle(" h2=1%, n=100, naffect=5")


load("D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/simulation/Simulation_script/script/ROC_sim_CpG_discovery/h2_01_n100_CpG_10.RData")
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
roc1 <- simple_roc(res$CpG, res$hmm_lfsr)
roc2 <- simple_roc(res$CpG,res$pv_ti)
roc3 <- simple_roc(res$CpG,res$affected_TI)
df_plot = data.frame (TPR= c(roc0$TPR,
                             roc1$TPR,
                             roc2$TPR ),
                      FPR=  c(roc0$FPR,
                              roc1$FPR,
                              roc2$FPR ),
                      col=factor(c( rep( "pv", length(roc0$FPR )),
                                    rep( "hmm", length(roc1$FPR )),
                                    rep( "TI_pv", length(roc1$FPR )) )
                      )
)

P12= ggplot(df_plot, aes(x=TPR, y=FPR, col=col))+geom_line()+
  theme_cowplot()+xlab("FDR")+ylab("Power")+
  ggtitle(" h2=1%, n=100, naffect=10")






load("D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/simulation/Simulation_script/script/ROC_sim_CpG_discovery/h2_05_n100_CpG_10_decay.RData")
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
roc1 <- simple_roc(res$CpG, res$hmm_lfsr)
roc2 <- simple_roc(res$CpG,res$pv_ti)
roc3 <- simple_roc(res$CpG,res$affected_TI)
df_plot = data.frame (TPR= c(roc0$TPR,
                             roc1$TPR,
                             roc2$TPR ),
                      FPR=  c(roc0$FPR,
                              roc1$FPR,
                              roc2$FPR  ),
                      col=factor(c( rep( "pv", length(roc0$FPR )),
                                    rep( "hmm", length(roc1$FPR )),
                                    rep( "TI_pv", length(roc1$FPR )) )
                      )
)

P21= ggplot(df_plot, aes(x=TPR, y=FPR, col=col))+geom_line()+
  theme_cowplot()+xlab("FDR")+ylab("Power")+
  ggtitle(" h2=5%, n=100, naffect=5")


load("D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/simulation/Simulation_script/script/ROC_sim_CpG_discovery/h2_05_n100_CpG_10.RData")
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
roc1 <- simple_roc(res$CpG, res$hmm_lfsr)
roc2 <- simple_roc(res$CpG,res$pv_ti)
roc3 <- simple_roc(res$CpG,res$affected_TI)
df_plot = data.frame (TPR= c(roc0$TPR,
                             roc1$TPR,
                             roc2$TPR ),
                      FPR=  c(roc0$FPR,
                              roc1$FPR,
                              roc2$FPR ),
                      col=factor(c( rep( "pv", length(roc0$FPR )),
                                    rep( "hmm", length(roc1$FPR )),
                                    rep( "TI_pv", length(roc1$FPR )) )
                      )
)

P22= ggplot(df_plot, aes(x=TPR, y=FPR, col=col))+geom_line()+
  theme_cowplot()+xlab("FDR")+ylab("Power")+
  ggtitle(" h2=5%, n=100, naffect=10")








load("D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/simulation/Simulation_script/script/ROC_sim_CpG_discovery/h2_10_n100_CpG_10_decay.RData")
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
roc1 <- simple_roc(res$CpG, res$hmm_lfsr)
roc2 <- simple_roc(res$CpG,res$pv_ti)
roc3 <- simple_roc(res$CpG,res$affected_TI)
df_plot = data.frame (TPR= c(roc0$TPR,
                             roc1$TPR,
                             roc2$TPR ),
                      FPR=  c(roc0$FPR,
                              roc1$FPR,
                              roc2$FPR ),
                      col=factor(c( rep( "pv", length(roc0$FPR )),
                                    rep( "hmm", length(roc1$FPR )),
                                    rep( "TI_pv", length(roc1$FPR )) )
                      )
)

P31= ggplot(df_plot, aes(x=TPR, y=FPR, col=col))+geom_line()+
  theme_cowplot()+xlab("FDR")+ylab("Power")+
  ggtitle(" h2=10%, n=100, naffect=5")


load("D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/simulation/Simulation_script/script/ROC_sim_CpG_discovery/h2_10_n100_CpG_10.RData")
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
roc1 <- simple_roc(res$CpG, res$hmm_lfsr)
roc2 <- simple_roc(res$CpG,res$pv_ti)
roc3 <- simple_roc(res$CpG,res$affected_TI)
df_plot = data.frame (TPR= c(roc0$TPR,
                             roc1$TPR,
                             roc2$TPR ),
                      FPR=  c(roc0$FPR,
                              roc1$FPR,
                              roc2$FPR ),
                      col=factor(c( rep( "pv", length(roc0$FPR )),
                                    rep( "hmm", length(roc1$FPR )),
                                    rep( "TI_pv", length(roc1$FPR )) )
                      )
)

P32= ggplot(df_plot, aes(x=TPR, y=FPR, col=col))+geom_line()+
  theme_cowplot()+xlab("FDR")+ylab("Power")+
  ggtitle(" h2=10%, n=100, naffect=10")



library(ggplot2)
library(gridExtra)
library(grid)
library(ggpubr)
library(cowplot)
library(dplyr)

# Create text labels for the column titles
titles <- list(
  textGrob(label = "WGBS = block", gp = gpar(fontsize = 16, fontface = "bold")),
  textGrob(label = "WGBS = distance decay", gp = gpar(fontsize = 16, fontface = "bold"))
)

# Create text labels for the row annotations (h^2 values) using LaTeX-style expressions
h2_labels <- list(
  textGrob(expression(h^2~"= 1%"), rot = 90, gp = gpar(fontsize = 14, fontface = "bold")),
  textGrob(expression(h^2~"= 5%"), rot = 90, gp = gpar(fontsize = 14, fontface = "bold")),
  textGrob(expression(h^2~"= 10%"), rot = 90, gp = gpar(fontsize = 14, fontface = "bold"))
)

# Extract the legend from one of the plots
legend_plot <- get_legend(P11 + theme(legend.position = "bottom", legend.justification = "center"))

# Arrange the grid with correct alignment
P_out <- grid.arrange(
  arrangeGrob(
    textGrob(""), titles[[1]], titles[[2]],  # Empty space for alignment
    ncol = 3, widths = c(0.08, 1, 1)
  ),
  arrangeGrob(
    h2_labels[[1]], P12 + theme(legend.position = "none"), P11 + theme(legend.position = "none"),
    ncol = 3, widths = c(0.08, 1, 1)
  ),
  arrangeGrob(
    h2_labels[[2]], P22 + theme(legend.position = "none"), P21 + theme(legend.position = "none"),
    ncol = 3, widths = c(0.08, 1, 1)
  ),
  arrangeGrob(
    h2_labels[[3]], P32 + theme(legend.position = "none"), P31 + theme(legend.position = "none"),
    ncol = 3, widths = c(0.08, 1, 1)
  ),
  arrangeGrob(legend_plot),
  heights = c(0.08, 1, 1, 1, 0.1)
)

# Display the final plot
P_out
