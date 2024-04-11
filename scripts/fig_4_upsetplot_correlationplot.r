library("dplyr")
library("readr")
library("stringr")
library("purrr")
library("tidyr")
library("ggplot2")
library("cowplot")
library("ComplexUpset")
fig1 = upset(read_delim("../data/fig_4_upsetplot_data//2_upset.tsv"),intersect = c("mQTL","haQTL","eQTL","pQTL"),
  keep_empty_groups = F,
      base_annotations=list(`Intersection size` = intersection_size( bar_number_threshold = 1, position = position_dodge(0.5), width = 0.3 ,text = list(size = 10)  )+ylab("")+annotate("point", y = 2500, x = c(4:14), color = "blue", size = 5  ) + annotate("point", y = 2500, x = c(7,8,13,14), color = "red", size = 5  )       ) ,
          themes=upset_default_themes(plot.margin=unit(c(0,0,0,20),"mm"),axis.text=element_text(size=40), axis.title.x = element_blank() ,text=element_text(size=40) )     ,  
          min_degree = 1)
fig2 = upset(read_delim("../data/fig_4_upsetplot_data//2_upset_superfine.tsv"),intersect = c("mQTL","haQTL","eQTL","pQTL"),
  keep_empty_groups = F,
      base_annotations=list(`Intersection size` = intersection_size( bar_number_threshold = 1, position = position_dodge(0.5), width = 0.3 ,text = list(size = 10)  )+ylab("")+annotate("point", y = 2500, x = c(4:14), color = "blue", size = 5  ) + annotate("point", y = 2500, x = c(6,8,13,14), color = "red", size = 5  )       ) ,
          themes=upset_default_themes(plot.margin=unit(c(0,0,0,20),"mm"),axis.text=element_text(size=40), axis.title.x = element_blank() ,text=element_text(size=40) )     ,  
          min_degree = 1)
fig1%>%ggsave(filename = "../plot//fig4_upsetplot.pdf", device = "pdf", width = 20, height = 15)
fig1%>%ggsave(filename = "../plot//fig4_upsetplot.png", device = "png", width = 20, height = 15)
fig2%>%ggsave(filename = "../plot//fig4_upsetplot_superfine.pdf", device = "pdf", width = 20, height = 15)
fig2%>%ggsave(filename = "../plot//fig4_upsetplot_superfine.png", device = "png", width = 20, height = 15)
cs_prop_16e_cat= read_delim("../data/fig_4_upsetplot_data//2_correlation.tsv","\t")
fig3= cs_prop_16e_cat%>%ggplot(aes(x = max_effect.ha, y = max_effect.m  ))+
    geom_smooth(aes(color = cat), size = 10,method = "lm")+geom_point(aes(color = cat))+
    geom_vline( xintercept = 0) + geom_hline(yintercept = 0) + 
    geom_abline(slope = -1) +xlab("Max haQTL effect")+ 
    ylab("Max mQTL effect")  + theme_bw()+   theme(text = element_text(size = 40))  +   scale_color_manual("Degree of affected TSS sharing",values = c("red","black","blue"))
cs_prop_16e_cat= read_delim("../data/fig_4_upsetplot_data//2_correlation_superfine.tsv","\t")
fig4= cs_prop_16e_cat%>%ggplot(aes(x = max_effect.ha, y = max_effect.m  ))+
    geom_smooth(aes(color = cat), size = 10,method = "lm")+geom_point(aes(color = cat))+
    geom_vline( xintercept = 0) + geom_hline(yintercept = 0) + 
    geom_abline(slope = -1) +xlab("Max haQTL effect")+ 
    ylab("Max mQTL effect")  + theme_bw()+   theme(text = element_text(size = 40))  +   scale_color_manual("Degree of affected TSS sharing",values = c("red","black","blue")) 
fig3%>%ggsave(filename = "../plot//fig4_correlation.pdf", device = "pdf", width = 20, height = 15)
fig3%>%ggsave(filename = "../plot//fig4_correlation.png", device = "png", width = 20, height = 15)
fig4%>%ggsave(filename = "../plot//fig4_correlation_superfine.pdf", device = "pdf", width = 20, height = 15)
fig4%>%ggsave(filename = "../plot//fig4_correlation_superfine.png", device = "png", width = 20, height = 15)