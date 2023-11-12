annotation_plot <- read.delim("D:/Document/Serieux/Travail/Data_analysis_and_papers/susiF_data_simu/enrichment_results/annotation_plot_df_firth_log_ha_m_context_grouped.tsv")
library(ggplot2)
library(dplyr)




t_pv <- 2*( 1-pnorm(abs(annotation_plot$mean), mean=0, sd=((annotation_plot$upper-annotation_plot$lower)/(1.96*2)) ) )
annotation_plot$pvalue  <- t_pv

sig <- rep( "Yes", length(annotation_plot$pvalue))
sig[which(annotation_plot$pvalue> 0.05/length(annotation_plot$pvalue))] <- "No"

annotation_plot$Significant <- factor( sig, levels = c("Yes","No"))

annotation_plot$feature<- sub("_"," ", annotation_plot$feature)
annotation_plot$feature<- sub("_"," ", annotation_plot$feature)

annotation_plot$feature<- sub("_"," ", annotation_plot$feature)
annotation_plot$feature<- sub("_"," ", annotation_plot$feature)

annotation_plot$feature<- sub("_"," ", annotation_plot$feature)
annotation_plot$feature<- sub("_"," ", annotation_plot$feature)
annotation_plot$feature<- sub("[.]"," ", annotation_plot$feature)
annotation_plot$feature<- sub("[.]"," ", annotation_plot$feature)

table(sig, annotation_plot$type)


ggplot(annotation_plot)+
  geom_point(aes(x = mean, y =  reorder(feature,-abs(mean) ), shape=Significant , color = type))+
  geom_linerange(aes(xmin = lower, xmax= upper , y = feature   ,color = type))+
  geom_vline(aes( xintercept = 0 ))+
  theme_bw()+
  ylab("")+
  theme(text = element_text(size = 15),
        axis.text.y= element_text(size = 7))+
  xlab("Log OR")+
  facet_grid(.~ind,scales = "free")+
  labs(color=NULL)

P1 <- ggplot(annotation_plot[ which(annotation_plot$ind=="conservation"),])+
  geom_point(aes(x = mean, y =  reorder(feature,-abs(mean)), shape=Significant, color = type))+
  geom_linerange(aes(xmin = lower, xmax= upper , y = feature ,  color = type))+
  geom_vline(aes( xintercept = 0 ))+
  theme_bw()+
  ylab("")+
  theme(text = element_text(size = 20) )+
  ggtitle("Conservation annotation")+
  xlab("Log OR")+
  labs(color=NULL)
P1
P2 <- ggplot(annotation_plot[ which(annotation_plot$ind=="dna_methylation"),])+
  geom_point(aes(x = mean, y =  reorder(feature,-abs(mean)), shape=Significant, color = type))+
  geom_linerange(aes(xmin = lower, xmax= upper , y = feature ,  color = type))+
  geom_vline(aes( xintercept = 0 ))+
  theme_bw()+
  ylab("")+
  theme(text = element_text(size = 20))+ggtitle("DNA methylation")+
  xlab("Log OR")+
  labs(color=NULL)
P2
P3 <- ggplot(annotation_plot[ which(annotation_plot$ind=="enhancer"),])+
  geom_point(aes(x = mean, y =  reorder(feature,-abs(mean)), shape=Significant, color = type))+
  geom_linerange(aes(xmin = lower, xmax= upper , y = feature ,  color = type))+
  geom_vline(aes( xintercept = 0 ))+
  theme_bw()+
  ylab("")+
  theme(text = element_text(size = 20) )+ggtitle("Enhancer  annotation")+
  xlab("Log OR")+
  labs(color=NULL)
P3

P4 <- ggplot(annotation_plot[ which(annotation_plot$ind=="genomic_features"),])+
  geom_point(aes(x = mean, y =  reorder(feature,-abs(mean)), shape=Significant, color = type))+
  geom_linerange(aes(xmin = lower, xmax= upper , y = feature ,  color = type))+
  geom_vline(aes( xintercept = 0 ))+
  theme_bw()+
  ylab("")+
  theme(text = element_text(size = 20))+
  xlab("Log OR")+
  labs(color=NULL)
P4



P5 <- ggplot(annotation_plot[ which(annotation_plot$ind=="histone_modifications"),])+
  geom_point(aes(x = mean, y =  reorder(feature,-abs(mean)), shape=Significant, color = type))+
  geom_linerange(aes(xmin = lower, xmax= upper , y = feature ,  color = type))+
  geom_vline(aes( xintercept = 0 ))+
  theme_bw()+
  ylab("")+
  theme(text = element_text(size = 20))+ggtitle("Histone annotation")+
  xlab("Log OR")+
  labs(color=NULL)
P5

P6 <- ggplot(annotation_plot[ which(annotation_plot$ind=="others"),])+
  geom_point(aes(x = mean, y =  reorder(feature,-abs(mean)), shape=Significant, color = type))+
  geom_linerange(aes(xmin = lower, xmax= upper , y = feature ,  color = type))+
  geom_vline(aes( xintercept = 0 ))+
  theme_bw()+
  ylab("")+
  theme(text = element_text(size = 20))+
  xlab("Log OR")+
  labs(color=NULL)
P6

P7 <- ggplot(annotation_plot[ which(annotation_plot$ind=="promoter"),])+
  geom_point(aes(x = mean, y =  reorder(feature,-abs(mean)), shape=Significant, color = type))+
  geom_linerange(aes(xmin = lower, xmax= upper , y = feature ,  color = type))+
  geom_vline(aes( xintercept = 0 )
  )+
  theme_bw()+
  ylab("")+
  theme(text = element_text(size = 20))+ggtitle("Promoter annotation")+
  xlab("Log OR")+
  labs(color=NULL)
P7



library(ggpubr)
ggarrange( P1,P3, P7,P5,
           common.legend = TRUE, legend="bottom")
