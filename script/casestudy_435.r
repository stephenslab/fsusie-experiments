library("dplyr")
library("readr")
library("stringr")
library("purrr")
library("tidyr")
library("ggplot2")
library("cowplot")

data_list = readRDS("../data/Case435.rds")
effect = data_list$effect
refine_plot = data_list$refine_plot
gene_plot = data_list$gene_plot 
annotation = data_list$annotation

n = c(1,2,3,4)
color  = c("black", "dodgerblue2", "#6A3D9A","#FF7F00","skyblue2","#6A3D9A",
                   "gold1",  "#FB9A99", "palegreen2",
                   "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1",
                   "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1",
                   "yellow4", "yellow3","darkorange4","brown","navyblue","#FF0000",
                   "darkgreen","#FFFF00","purple","#00FF00","pink","#0000FF",
                   "orange","#FF00FF","cyan","#00FFFF","#FFFFFF")

### Effect plot

refine_effect_plot_plot<-  ggplot( effect  )+
  geom_line(aes(x = pos, y = mid,color = effect), size=1.4)+
  geom_ribbon(aes(x = pos,
                  y = mid,
                  ymin = low,
                  ymax = up,
                  fill = effect),
              alpha = 0.2)+
  ylab("Estimated Effect") +
  xlab("POS")+
  facet_grid(factor(type,levels =  c("haQTL_Effect_1",
                                     "mQTL_Effect_1",
                                     "haQTL_Effect_3",
                                     "mQTL_Effect_3"),
                    labels =   c("haQTL CS 1",
                                 "mQTL CS 1",
                                 "haQTL CS 3",
                                 "mQTL CS 3")
  )~., scales = "free" )+
  scale_color_manual("Credible set",values = color[n+1])+
  geom_hline(aes(yintercept = 0), color = "black")+
  scale_fill_manual("Credible set",values = color[n+1])+
  theme_bw()+
  xlab("") +
  ylab("Estimated Effect")+
  theme(text = element_text(size = 10),
        legend.position="none")+
  theme(plot.margin=unit(c(5,0,0,0),"mm"),
        strip.text.y.right = element_text(angle = 0,
                                          size = 10),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10)
  )
refine_effect_plot_plot

### PIP plot

color2 = c("black", "dodgerblue2","#FF7F00", "#6A3D9A","skyblue2",
                   "gold1",  "#FB9A99", "palegreen2",
                   "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1",
                   "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1",
                   "yellow4", "yellow3","darkorange4","brown","navyblue","#FF0000",
                   "darkgreen","#FFFF00","purple","#00FF00","pink","#0000FF",
                   "orange","#FF00FF","cyan","#00FFFF","#FFFFFF")

refine_plot_plot  <-  ggplot2::ggplot(refine_plot,aes(y = y,
                                                 x = pos,
                                                 col =  as.factor(new_CS),
                                                 shape = Shared )) +
    geom_text( aes(x =pos + 6000, y = y ,color =  as.factor(new_CS), label = sign ),alpha = 0.5 , size = 10)+
  facet_grid(molecular_trait_id ~.)+
  geom_point(size = 4) +
  scale_color_manual("CS",values = color2) +
  theme_bw()+
  theme(axis.ticks.x = element_blank()) +
  ylab("Posterior Inclusion Probability (PIP)")+
  xlim(plot_range)+
  theme(axis.ticks.x = element_blank() ) +
  theme(strip.text.y.right = element_text(angle = 0))+
  xlab("") +
  theme(text = element_text(size = 10),legend.position = "bottom")


### Gene plot

plot_range = c(153120464,153679654)
gene_plot$x_label <- (0.5*(gene_plot$end-gene_plot$start)+gene_plot$start)
#gene_plot$x_label[which(gene_plot$gene_name=="CHRM5")] <- gene_plot$x_label[which(gene_plot$gene_name=="CHRM5")]+20000
n = 0.9
gene_plot_plot <-ggplot(gene_plot,aes()) +
  geom_segment( aes(x = start,xend = end, y = 0.88, yend = 0.88 ) ,
                arrow = arrow(length = unit(0.5, "cm")) )+
  geom_text(aes(x = x_label,y = 0.885, label = gene_name))+
  ylab("")+
  xlab("")+
  theme(legend.position="none")+
  theme(text = element_blank())+
  xlim(plot_range)+theme_bw()+
                                     theme(axis.text.x = element_blank(),
                                           axis.text.y = element_blank()) +
                                     theme(strip.text.y.right = element_text(angle = 0))+
                                     xlab("") +geom_point( aes(x = X4,xend = X5, y = 0.88 , yend = 0.88 ),size = 2, color  = "#6A3D9A", data = annotated%>%count(X3,X4,X5)%>%filter(X3 %in% c("transcript")))+
                                     xlab("") +geom_point( aes(x = X4,xend = X5, y = 0.88 , yend = 0.88 ),size = 2, color  = "blue", data = annotated%>%count(X3,X4,X5)%>%filter(X3 %in% c("transcript"), X4 ==153601955 |X4 ==153632364  ))+
                                   geom_segment( aes(x = start,xend = end, y = 0.87, yend = 0.87 ),alpha = 0.7,size = 5, color  = "#6A3D9A", data = annotation%>%select(start = start1, end = end1, chr = chr1)%>%filter(chr == "chr4", start %in% c(153333848,153343848)  )  ) +

                                   geom_segment( aes(x = start,xend = end, y = 0.87, yend = 0.87 ),alpha = 0.7,size = 5, color  = "#6A3D9A", data = annotation%>%select(start = start2, end = end2, chr = chr2,start1)%>%filter(chr == "chr4", start1  %in% c(153333848,153343848) )  ) +
                                    geom_segment( aes(x = start1,xend = (start1+start2)/2, y = 0.87, yend = 0.835 ),alpha = 0.6,size = 0.5, color  = "#6A3D9A", data = annotation%>%filter(chr1 == "chr4", start1  == 153333848    )  ) +
                                   geom_segment( aes(x = (start1+start2)/2,xend = start2, y = 0.835, yend = 0.87 ),alpha = 0.6,size = 0.5, color  = "#6A3D9A", data = annotation%>%filter(chr1 == "chr4", start1  == 153333848   )  ) +
                                   geom_segment( aes(x = start1,xend = (start1+start2)/2, y = 0.87, yend = 0.835 ),alpha = 0.6,size = 0.5, color  = "#6A3D9A", data = annotation%>%filter(chr1 == "chr4", start1  == 153343848   )  ) +
                                   geom_segment( aes(x = (start1+start2)/2,xend = start2, y = 0.835, yend = 0.87 ),alpha = 0.6,size = 0.5, color  = "#6A3D9A", data = annotation%>%filter(chr1 == "chr4", start1  == 153343848  )  ) +
                                  
                                   xlab("") 

### Combined:

cowplot::plot_grid(plotlist = list(refine_effect_plot_plot+theme(strip.text.y.right = element_text(angle = 0,size = 10),panel.spacing=unit(0.7, "lines"),axis.text.y = element_text(size = 10))+
                                     xlim(plot_range)+
                                     theme(text = element_text(size = 10)),#microglia_enhancer_activity,#+facet_grid(TargetGene_name~., scales = "free"),
                                    gene_plot_plot +
                                   
                                     theme(text = element_text(size = 10)),
                                   refine_plot_plot+
                                     theme_bw()+
                                     theme(axis.ticks.x = element_blank()) +
                                     theme(strip.text.y.right = element_text(angle = 0,size = 10))+
                                     xlab("") +ylim(c(0,1))+
                                     theme(text = element_text(size = 10),axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), 
                                           panel.spacing=unit(0.7, "lines"),legend.position = "none")+scale_y_continuous(breaks = c(0,0.5,1))
                                   
) ,
ncol = 1,
align = "v",
axis = "tlbr",
rel_heights = c(5,3,6),labels  = c("A","B","C"),label_size = 10
) -> result_plot
result_plot
result_plot%>%ggsave(filename = "../plot/casestudy_435_highlight.pdf", device = "pdf",
       width = 29.7,
       height = 21,
       units = "cm"
)
result_plot%>%ggsave(filename = "../plot/casestudy_435_highlight.png",device = "png",
       width = 29.7,
       height = 21,
       units = "cm"
)