rm(list=ls())
library("dplyr")
library("readr")
library("stringr")
library("purrr")
library("tidyr")
library("ggplot2")
library("cowplot")
library("ComplexUpset")
path <- getwd()
data_list =  readRDS(paste(path, "/data/Case1411.rds", sep=""))
effect = data_list$effect
refine_plot = data_list$refine_plot
gene_plot = data_list$gene_plot 
annotation = data_list$annotation
tf_candidates = data_list$tf_candidates
plot_range = c(35930000,36530000)
color2 = c("black", "dodgerblue2","#FF7F00", "#6A3D9A","skyblue2",
                   "gold1",  "#FB9A99", "palegreen2",
                   "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1",
                   "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1",
                   "yellow4", "yellow3","darkorange4","brown","navyblue","#FF0000",
                   "darkgreen","#FFFF00","purple","#00FF00","pink","#0000FF",
                   "orange","#FF00FF","cyan","#00FFFF","#FFFFFF")

### Effect plot
n = c(1,2,3,4)
nn = 0.9
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
                                     "mQTL_Effect_1"
                                     ),
                    labels =   c("haQTL CS 1",
                                 "mQTL CS 1")
  )~., scales = "free" )+
  scale_color_manual("Credible set",values = color2[n+1])+
  geom_hline(aes(yintercept = 0), color = "black")+
  scale_fill_manual("Credible set",values = color2[n+1])+
  theme_bw()+
  xlab("") +
  ylab("")+
  theme( 
    legend.position="none",
    plot.margin=unit(c(5,0,0,0),"mm"),
    strip.text.y.right = element_text(angle = 0,
                                      size = 10),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),
    panel.border = element_rect(colour = "black", fill=NA, size=1.2),
    strip.background =element_rect(fill="white"),
    
    panel.spacing=unit(0.7, "lines")
  )+ xlim(plot_range) 
 

refine_effect_plot_plot

### PIP plot

refine_plot$molecular_trait_id <- factor(refine_plot$molecular_trait_id , levels= c("DNAm",
                                                                                    "H3K9ac",
                                                                                    "eQTL APOL2",
                                                                                    "pQTL APOL2",
                                                                                    "Inh APOL2",
                                                                                    "Neu eQTL APOL2" 
))

refine_plot_plot  <-  ggplot2::ggplot(refine_plot,aes(y = y,
                                                 x = pos,
                                                 col =  as.factor(new_CS),
                                                 shape = Shared )) +
  geom_text( aes(x =pos + 6000, y = y ,color =  as.factor(new_CS), label = sign ),alpha = 0.5 , size = 10)+
  facet_grid(molecular_trait_id ~.)+
  geom_point(size = 4) +
  scale_color_manual("CS",values = color2) +
  theme_bw()+ 
  ylab("")+
  xlim(plot_range)+ 
  theme_bw()+
  theme(axis.ticks.x = element_blank() ,
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10), 
        panel.spacing=unit(0.7, "lines") , 
        panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        strip.background =element_rect(fill="white"),
        legend.position = "none",
        strip.text.y.right = element_text(angle = 0,size = 10),
        text = element_text(size = 10) )+
  xlab("") +ylim(c(0,1))+  
  scale_y_continuous(breaks = c(0, 1), limits =c(0,1.25) )   
refine_plot_plot
### TF plot

tf_plot = ggplot(data = annotation%>%arrange(start1)%>%filter(start2 > 36500000))+
  geom_segment(aes(color = region ,
                   x = start1,
                   xend = end1,
                   y = 0.87,
                   yend =  0.87),
               alpha = 0.7,
               size = 5)+
  geom_segment(aes(x = start2,
                   xend = end2,
                   y = 0.87, 
                   yend =  0.87),
               alpha = 0.7,
               size = 5,
               color = "green")+
  geom_segment(aes( x = start1,
                    xend = (start1+start2)/2,
                    y = 0.87 , 
                    yend = 0.835 )
               )+
  geom_segment(aes( x = (start1+start2)/2 ,
                    xend = start2, 
                    y = 0.835 , 
                    yend = 0.87 )
               )+
  geom_segment(data = tf_candidates%>%
                      filter() ,
               aes(x = X2,
                  xend = X3,
                   y = 1.1-X7/30,
                  yend = 1.1-X7/30),# bars
               size = 2
               )+
  geom_segment( aes(x = start,
                    xend = end, 
                    y = (nn-strand/100),
                    yend =(nn-strand/100)
                    ) ,
                arrow = arrow(length = unit(0.5, "cm")),
                data = gene_plot # arrow gene
                )+
  geom_text(aes(x = 35985000, 
                y = 1.1-X7/30,
                label = X4 ),
            size = 5,check_overlap = TRUE,
            data =tf_candidates%>%filter() # names of the TF
            )+
  geom_text(data = gene_plot,
            aes(x = (x_label+30000),
                y = 0.885, 
                label = gene_name,
                vjust=-1),#name gene
            size = 5
            )+
  scale_color_manual("region",
                     values = c("orange",
                                "red",
                                "yellow",
                                "purple") 
                     )+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank() ,
        axis.ticks.y = element_blank() )+
  ylab("")+
  xlim(plot_range[1],plot_range[2])
tf_plot

### Combined:

cowplot::plot_grid(plotlist = list(refine_effect_plot_plot,
                                    tf_plot  ,
                                   refine_plot_plot 
                                   
) ,
ncol = 1,
align = "v",
axis = "tlbr",
rel_heights = c(4,4,6),labels  = c("A","B","C"),label_size = 10
) -> result_plot
result_plot

result_plot%>%ggsave(filename = paste(path, "/plot/casestudy_1411_highlight.pdf",sep=""), device = "pdf",
       width = 29.7,
       height = 21,
       units = "cm"
)
 