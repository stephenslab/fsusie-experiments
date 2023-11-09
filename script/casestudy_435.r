library("dplyr")
library("readr")
library("stringr")
library("purrr")
library("tidyr")
library("ggplot2")
library("cowplot")
library("ComplexUpset")
corresponding_vector <- c(153483848, 153488848, 153498848, 153523848, 153528848, 153533848, 153563848, 153573848, 153583848)

# Not Corresponding to 153343848 in Start1
not_corresponding_vector <- c(153488848, 153523848, 153528848, 153533848, 153573848)

n = c(1,2,3,4)
read_delim("../data/resource/Homo_sapiens.GRCh38.103.chr.reformatted.ERCC_chr4_chr22.gtf.gz", col_names   = 0, skip = 1)%>%
filter(X1 == "chr4", X4 > plot_range[[1]],X5 < plot_range[[2]] , str_detect(X9, "TRIM2" ) |str_detect(X9, "TMEM131L" ) | str_detect(X9, "MND1" )   )%>%
mutate(gene_name =  map(X9,~read.table(text = .x , sep = "\"")$V6 ))-> annotated

color = color2 = c("black", "dodgerblue2", "#6A3D9A","#FF7F00","skyblue2","#6A3D9A",
                   "gold1",  "#FB9A99", "palegreen2",
                   "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1",
                   "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1",
                   "yellow4", "yellow3","darkorange4","brown","navyblue","#FF0000",
                   "darkgreen","#FFFF00","purple","#00FF00","pink","#0000FF",
                   "orange","#FF00FF","cyan","#00FFFF","#FFFFFF")



effect <- read_delim("../data/case_study_tad435/4_435_effect_refined.tsv")

#refine_effect_plot =  read_delim(paste(path,"4_1182_effect.tsv", sep=""))
refine_effect_plot<-  ggplot( effect  )+
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
  scale_color_manual("Credible set",values = color2[n+1])+
  geom_hline(aes(yintercept = 0), color = "black")+
  scale_fill_manual("Credible set",values = color2[n+1])+
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


refine_effect_plot
gene_plot  = read_delim("../data/case_study_tad435/4_435_gene.tsv")
#gene_plot[4,] = gene_plot[3,]
#gene_plot[4,] = list("chr15",34075155,33858782, "ENSG00000169857" , "AVEN" ,2,2 )
plot_range = c(153120464,153679654)
gene_plot$x_label <- (0.5*(gene_plot$end-gene_plot$start)+gene_plot$start)
#gene_plot$x_label[which(gene_plot$gene_name=="CHRM5")] <- gene_plot$x_label[which(gene_plot$gene_name=="CHRM5")]+20000
n = 0.9
gene_plot_plot <-ggplot(gene_plot,aes()) +
  geom_segment( aes(x = start,xend = end, y = (n-strand/100), yend =(n-strand/100) ) ,
                arrow = arrow(length = unit(0.5, "cm")) )+
  geom_text(aes(x = x_label,y = (n-0.05-strand/100), label = gene_name, vjust=-1),
            size = 7)+
  ylab("")+
  xlab("")+
  theme(legend.position="none")+
  theme(text = element_blank())+
  xlim(plot_range)+
  geom_point(aes(x = start, y = (n-strand/100)),
             color = "black",size = 3  ) 


refine_plot =read_delim("../data/case_study_tad435/4_435_pip_oli.tsv")
refine_plot$molecular_trait_id  [ which(refine_plot$molecular_trait_id =="Oli_eQTL_TMEM131L" )]<- "Oli TMEM131L"
refine_plot$molecular_trait_id  [ which(refine_plot$molecular_trait_id =="Oli_eQTL_MND1" )]<- "Oli MND1"
refine_plot$molecular_trait_id  [ which(refine_plot$molecular_trait_id =="Oli_eQTL_TRIM2" )]<- "Oli TRIM2"
refine_plot$molecular_trait_id  [ which(refine_plot$molecular_trait_id =="eQTL_TMEM131L" )]<- "TMEM131L"
refine_plot$molecular_trait_id  [ which(refine_plot$molecular_trait_id =="haQTL_tad435" )]<- "H3K9ac"
refine_plot$molecular_trait_id  [ which(refine_plot$molecular_trait_id =="mQTL_tad435" )]<- "DNAm"
refine_plot$molecular_trait_id <- factor(refine_plot$molecular_trait_id , levels = c("Oli MND1","Oli TRIM2", "Oli TMEM131L","TMEM131L",
                                                                                      "H3K9ac",
                                                                                      "DNAm" ))
color = color2 = c("black", "dodgerblue2","green4","#6A3D9A","skyblue2", "#FF7F00",
                   "gold1",  "#FB9A99", "palegreen2",
                   "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1",
                   "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1",
                   "yellow4", "yellow3","darkorange4","brown","navyblue","#FF0000",
                   "darkgreen","#FFFF00","purple","#00FF00","pink","#0000FF",
                   "orange","#FF00FF","cyan","#00FFFF","#FFFFFF")%>%unique()
#### Change coloring of non-overlap CS
refine_plot%>%mutate(new_CS = case_when(molecular_trait_id == "Oli TMEM131L" & new_CS == 2 ~  9,
                                    molecular_trait_id == "H3K9ac" & new_CS == 2 ~  10,
                                     molecular_trait_id == "DNAm" & !Shared & new_CS == 1 ~  11,
                                    .default = new_CS 
                                   ))-> refine_plot
refine_plot%>%mutate(new_CS = case_when(
                                     molecular_trait_id == "DNAm" & !Shared & CS == 3 ~  12,
                                    .default = new_CS 
                                   ))-> refine_plot



### Add sign annotation, see pseudobulk notebook for code to gather this info
sign = 
 refine_plot%>%group_by(molecular_trait_id, new_CS)%>%
    summarize(pos = max(pos[which(y == max(y))]), y = max(y))%>%filter(new_CS != 0)%>%
    filter(!molecular_trait_id%in% c("H3K9ac", "DNAm")) %>%ungroup%>%mutate(sign = c("-","-","-","+","-") )
refine_plot = left_join(refine_plot,sign)



table(refine_plot$molecular_trait_id )
refine_plot  <-  ggplot2::ggplot(refine_plot,aes(y = y,
                                                 x = pos,
                                                 col =  as.factor(new_CS),
                                                 shape = Shared )) +
    geom_text( aes(x =pos + 6000, y = y ,color =  as.factor(new_CS), label = sign ),alpha = 0.5 , size = 10)+
  facet_grid(molecular_trait_id ~.)+
  geom_point(size = 4) +
  scale_color_manual("CS",values = color) +
  theme_bw()+
  theme(axis.ticks.x = element_blank()) +
  ylab("Posterior Inclusion Probability (PIP)")+
  xlim(plot_range)+
  theme(axis.ticks.x = element_blank() ) +
  theme(strip.text.y.right = element_text(angle = 0))+
  xlab("") +
  theme(text = element_text(size = 10),legend.position = "bottom")

cowplot::plot_grid(plotlist = list(refine_effect_plot+theme(strip.text.y.right = element_text(angle = 0,size = 20),panel.spacing=unit(0.7, "lines"),axis.text.y = element_text(size = 18))+
                                     xlim(plot_range)+
                                     theme(text = element_text(size = 30)),#microglia_enhancer_activity,#+facet_grid(TargetGene_name~., scales = "free"),
                                   gene_plot_plot+
                                     theme_bw()+
                                     theme(axis.text.x = element_blank(),
                                           axis.text.y = element_blank()) +
                                     theme(strip.text.y.right = element_text(angle = 0))+
                                     xlab("") +geom_point( aes(x = X4,xend = X5, y = 0.89 , yend = 0.89 ),size = 2, color  = "darkgreen", data = annotated%>%count(X3,X4,X5)%>%filter(X3 %in% c("transcript")))+
                                     xlab("") +geom_point( aes(x = X4,xend = X5, y = 0.89 , yend = 0.89 ),size = 2, color  = "red", data = annotated%>%count(X3,X4,X5)%>%filter(X3 %in% c("transcript"), X4 ==153601955 |X4 ==153632364  ))+
                                   #geom_segment( aes(x = X4,xend = X5, y = 0.815 , yend = 0.815 ),size = 7, color  = "red", data = annotated%>%count(X3,X4,X5)%>%filter(X3 %in% c("exom")))+
                                     #geom_point( aes(x = start,xend = end, y = 0.815 , yend = 0.815 ),size = 1, color  = "blue", data = cpg%>%filter(`#chr` == "chr4"))+
                                    # geom_point( aes(x = 153297198 ,xend = end, y = 0.815 , yend = 0.815 ),size = 1, color  = "red", data = cpg%>%filter(`#chr` == "chr4"))+
                                   ##geom_point( aes(x = start,xend = end, y = 0.815 , yend = 0.815),alpha = 0.5,size = 2, color  = "green", data = ha%>%filter(`#chr` == "chr4")) +
                                     #geom_segment( aes(x = start,xend = end, y = 0.5, yend = 0.5 ),size = 10, color  = "yellow", data = annot_ATAT_seq%>%filter(chr == "chr15", start > plot_range[[1]],end < plot_range[[2]]  )  ) +
                                     #geom_segment( aes(x = start,xend = end, y = 0.825, yend = 0.825 ),size = 3,color  = "cyan", data = read_delim("./added_annotation/Oli_hg38.pro.Nott.tsv")%>%filter(chr == "chr15", start > plot_range[[1]],end < plot_range[[2]]  )  ) +
                                   geom_segment( aes(x = start,xend = end, y = 0.88, yend = 0.88 ),alpha = 0.7,size = 5, color  = "cyan", data = read_delim("../data/resource//Oli_hg38.Int.Nott.tsv")%>%select(start = start1, end = end1, chr = chr1)%>%filter(chr == "chr4", start > plot_range[[1]],end < plot_range[[2]]  )  ) +
                                   geom_segment( aes(x = start,xend = end, y = 0.835, yend = 0.835 ),alpha = 0.7,size = 5, color  = "cyan", data = read_delim("../data/resource//Oli_hg38.Int.Nott.tsv")%>%select(start = start2, end = end2, chr = chr2,start1)%>%filter(chr == "chr4", start1 > plot_range[[1]],end < plot_range[[2]]  )  ) +
                                   #geom_segment( aes(x = start,xend = end, y = 0.835, yend = 0.835 ),alpha = 0.7,size = 5, color  = "blue", data = read_delim("./added_annotation/Oli_hg38.Int.Nott.tsv")%>%select(start = start1, end = end1, chr = chr1)%>%filter(chr == "chr4", start %in% c(153333848,153343848)  )  ) +

                                   #geom_segment( aes(x = start,xend = end, y = 0.88, yend = 0.88 ),alpha = 0.7,size = 5, color  = "blue", data = read_delim("./added_annotation/Oli_hg38.Int.Nott.tsv")%>%select(start = start2, end = end2, chr = chr2,start1)%>%filter(chr == "chr4", start1  %in% c(153333848,153343848) )  ) +
                                   #geom_segment( aes(x = start,xend = end, y = 0.835, yend = 0.835 ),alpha = 0.7,size = 5, color  = "cyan", data = read_delim("./added_annotation/Ast_hg38.Enh.Nott.tsv")%>%filter(chr == "chr4", start > plot_range[[1]],end < plot_range[[2]]  )  ) +
                                   #geom_segment( aes(x = start,xend = end, y = 0.835, yend = 0.835 ),alpha = 0.7,size = 5, color  = "blue", data = read_delim("./added_annotation/microglia_hg38.Enh.Nott.tsv")%>%filter(chr == "chr4", start > plot_range[[1]],end < plot_range[[2]]  )  ) +
                                   #geom_segment( aes(x = start,xend = end, y = 0.835, yend = 0.835 ),alpha = 0.7,size = 5, color  = "red", data = read_delim("./added_annotation/Neu_hg38.Enh.Nott.tsv")%>%filter(chr == "chr4", start > plot_range[[1]],end < plot_range[[2]]  )  ) +
                                   geom_segment( aes(x = start1,xend = start2, y = 0.88, yend =  0.835),alpha = 0.2,size = 0.5, color  = "cyan", data = read_delim("../data/resource/Oli_hg38.Int.Nott.tsv")%>%filter(chr1 == "chr4", start2 > plot_range[[1]],end2 < plot_range[[2]]  )  ) +
                                   geom_segment( aes(x = start1,xend = start2, y = 0.88, yend = 0.835 ),alpha = 0.6,size = 0.5, color  = "purple", data = read_delim("../data/resource/Oli_hg38.Int.Nott.tsv")%>%filter(chr1 == "chr4", start1  == 153333848 , start2  %in% corresponding_vector   )  ) +
                                   geom_segment( aes(x = start1,xend = start2, y = 0.88, yend = 0.835 ),alpha = 0.6,size = 0.5, color  = "purple", data = read_delim("../data/resource/Oli_hg38.Int.Nott.tsv")%>%filter(chr1 == "chr4", start1  == 153343848 , start2  %in% not_corresponding_vector   )  ) +
                                   xlab("") +geom_segment( aes(x = X4,xend = X5, y = 0.89 , yend = 0.89 ),size = 4,alpha = 0.3, color  = "green", data = read_delim("../data/case_study_tad435/4_435_gene_exom.tsv","\t"))+  
                                   #geom_segment(aes(x = 153234136,xend = 153248662, y = 0.9,yend = 0.9),color = "red")+ geom_segment(aes(x = 153295613,xend = 153307939, y = 0.9,yend = 0.9),color = "red") +
                                   #geom_point( aes(x = BP, y = 0.6 ),size = 1, color  = "black", data = enh)+
                                     theme(text = element_text(size = 20)),
                                   refine_plot+
                                     theme_bw()+
                                     theme(axis.ticks.x = element_blank()) +
                                     theme(strip.text.y.right = element_text(angle = 0,size = 20))+
                                     xlab("") +ylim(c(0,1))+
                    
                        
                                     theme(text = element_text(size = 20),axis.text.x = element_text(size = 30), axis.text.y = element_text(size = 18), 
                                           panel.spacing=unit(0.7, "lines"),legend.position = "none")+scale_y_continuous(breaks = c(0,0.5,1))
                                   
) ,
ncol = 1,
align = "v",
axis = "tlbr",
rel_heights = c(6,2,8),labels  = c("A","B","C"),label_size = 25
) -> result_plot
result_plot%>%ggsave(filename = "../plot/casestudy_435.pdf",width = 20, height = 15, device = "pdf")
result_plot%>%ggsave(filename = "../plot/casestudy_435.png",width = 20, height = 15, device = "png")
