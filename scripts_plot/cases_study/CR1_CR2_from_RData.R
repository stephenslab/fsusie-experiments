library(data.table)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(cowplot)
source("get_gene_annotations.R")
source("interpolate_effect_estimates.R")
load("../../outputs/CR1_CR2_obj.RData")
gene_file <-
  file.path("../../data/genome_annotations",
    "Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz")
genes <- get_gene_annotations(gene_file)

pos0 <- 207.45e6
pos1 <- 207.75e6
key_marker <- 207.577223

# The top panel shows the Alzheimer's Disease (AD) association
# p-values.
pdat1 <- as.data.frame(obj_plot$plot_df)
pdat1 <- subset(pdat1,
                study == "AD_Bellenguez_2022" &
                pos >= pos0 &
                pos <= pos1)
pdat1 <- pdat1[c("variant_alternate_id","pos","CS1","-log10(P)")]
pdat1 <- transform(pdat1,pos = pos/1e6)
names(pdat1) <- c("id","pos","CS","pval")
ids <- pdat1$id
ids[] <- NA
ids[pdat1$id == "chr1:207510847:T:G"] <- "rs12037841"
ids[pdat1$id == "chr1:207577223:T:C"] <- "rs679515"
ids[pdat1$id == "chr1:207623552:A:T"] <- "rs10863417"
ids[pdat1$id == "chr1:207624893:C:G"] <- "rs10863418"
ids[pdat1$id == "chr1:207629207:A:C"] <- "rs4844610"
pdat1$id <- ids
# > subset(pdat1,CS)
#                 id   pos   CS  pval
# chr1:207510847:T:G 207.5 TRUE 31.75
# chr1:207577223:T:C 207.6 TRUE 32.25
# chr1:207623552:A:T 207.6 TRUE 30.59
# chr1:207624893:C:G 207.6 TRUE 30.67
# chr1:207629207:A:C 207.6 TRUE 31.95
p1 <- ggplot(pdat1,aes(x = pos,y = pval,color = CS,label = id)) +
  geom_point(size = 0.75) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  geom_text_repel(size = 2.25,color = "dimgray",segment.color = "dimgray",
                  min.segment.length = 0,max.overlaps = Inf) +
  scale_color_manual(values = c("black","dodgerblue")) +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(207,208,0.05),
                     labels = NULL) +
  ylim(0,45) +
  labs(x = "",y = "AD") + 
  theme_cowplot(font_size = 9)

# The second panel shows the CR1 eQTL p-values.
pdat2 <- as.data.frame(obj_plot$plot_df)
pdat2 <- subset(pdat2,
                study == "DLPFC_DeJager_eQTL" &
                region == "CR1" &
                pos >= pos0 &
                pos <= pos1)
pdat2$CS <- as.numeric(NA)
pdat2[which(pdat2$CS1),"CS"] <- 1
pdat2[which(pdat2$CS3),"CS"] <- 3
pdat2 <- pdat2[c("variant_alternate_id","pos","CS","-log10(P)")]
names(pdat2) <- c("id","pos","CS","pval")
pdat2 <- transform(pdat2,
                   pos = pos/1e6,
                   CS  = factor(CS))
#
# > subset(pdat2,CS == 1)
#                 id   pos   CS  pval
# chr1:207564732:T:C 207.6 TRUE 64.56
# chr1:207573951:A:G 207.6 TRUE 63.47
# chr1:207577223:T:C 207.6 TRUE 64.89
# chr1:207611623:A:G 207.6 TRUE 64.47
# chr1:207612944:A:G 207.6 TRUE 64.47
# chr1:207613197:A:G 207.6 TRUE 64.51
# chr1:207613483:A:G 207.6 TRUE 64.47
# chr1:207629207:A:C 207.6 TRUE 62.42
#
# > nrow(subset(pdat2,CS == 3))
# 12
#
ids <- pdat2$id
ids[] <- NA
ids[pdat2$id == "chr1:207577223:T:C"] <- "rs679515"
pdat2$id <- ids
p2 <- ggplot(pdat2,aes(x = pos,y = pval,color = CS,label = id)) +
  geom_point(size = 0.75) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  scale_color_manual(values = c("limegreen","gold"),na.value = "black") +
  geom_text_repel(size = 2.25,color = "dimgray",segment.color = "dimgray",
                  min.segment.length = 0,max.overlaps = Inf) +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(207,208,0.05),
                     labels = NULL) +
  labs(x = "",y = "CR1 eQTL") + 
  theme_cowplot(font_size = 9)

# The third panel shows the CR2 eQTL p-values.
pdat3 <- as.data.frame(obj_plot$plot_df)
pdat3 <- subset(pdat3,
                study == "DLPFC_DeJager_eQTL" &
                region == "CR2" &
                pos >= pos0 &
                pos <= pos1)
pdat3 <- pdat3[c("variant_alternate_id","pos","CS1","-log10(P)")]
names(pdat3) <- c("id","pos","CS","pval")
pdat3 <- transform(pdat3,pos = pos/1e6)
#
# > subset(pdat2,CS == 1)
#                 id   pos   CS  pval
# chr1:207564732:T:C 207.6 TRUE 64.56
# chr1:207573951:A:G 207.6 TRUE 63.47
# chr1:207577223:T:C 207.6 TRUE 64.89
# chr1:207611623:A:G 207.6 TRUE 64.47
# chr1:207612944:A:G 207.6 TRUE 64.47
# chr1:207613197:A:G 207.6 TRUE 64.51
# chr1:207613483:A:G 207.6 TRUE 64.47
# chr1:207629207:A:C 207.6 TRUE 62.42
#
# > nrow(subset(pdat3,CS == 3))
# 12
#
ids <- pdat3$id
ids[] <- NA
# ids[pdat3$id == "chr1:207577223:T:C"] <- "rs679515"
pdat3$id <- ids
p3 <- ggplot(pdat3,aes(x = pos,y = pval,color = CS,label = id)) +
  geom_point(size = 0.75) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  scale_color_manual(values = c("black","darkorchid"),na.value = "black") +
  geom_text_repel(size = 2.25,color = "dimgray",segment.color = "dimgray",
                  min.segment.length = 0,max.overlaps = Inf) +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(207,208,0.05),
                     labels = NULL) +
  labs(x = "",y = "CR2 eQTL") + 
  theme_cowplot(font_size = 9)

# The fourth panel shows the haSNP PIPs.
ids   <- names(obj_plot$pip_fsusie_obj)
pdat4 <- data.frame(id  = as.character(NA),
                   pos = sapply(strsplit(ids,":",fixed = TRUE),"[[",2),
                   pip = obj_plot$pip_fsusie_obj,
                   cs  = as.character(NA),
                   stringsAsFactors = FALSE)
pdat4 <- transform(pdat4,pos = as.numeric(pos))
n <- length(obj_plot$cs_fsusie_obj)
for (i in 1:n) {
  snps <- names(obj_plot$cs_fsusie_obj[[i]])
  pdat4[snps,"cs"] <- i
}
pdat4 <- subset(pdat4,pos >= pos0 & pos <= pos1)
pdat4 <- transform(pdat4,pos = pos/1e6)
# > subset(pdat4,cs == 5 & pip > 0.05)
#                           id   pos    pip cs
# chr1:207577223:T:C      <NA> 207.6 0.1375  5
# chr1:207598421:CT:CTT   <NA> 207.6 0.3343  5
# chr1:207619376:CAAA:CAA <NA> 207.6 0.2049  5
pdat4[c("chr1:207577223:T:C",
        "chr1:207598421:CT:CTT",
        "chr1:207619376:CAAA:CAA"),"id"] <-
  c("rs679515","rs1168807665","rs869302047")
p4 <- ggplot(pdat4,aes(x = pos,y = pip,color = cs,label = id)) +
  geom_point(size = 0.75) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  geom_text_repel(size = 2.25,color = "dimgray",segment.color = "dimgray",
                  min.segment.length = 0,max.overlaps = Inf) +
  scale_color_manual(values = c("tomato","darkorange"),na.value = "black") +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(207,208,0.05),
                     labels = NULL) +
  ylim(0,0.4) +
  labs(x = "",y = "haSNP PIP") + 
  theme_cowplot(font_size = 9)

# The fifth panel shows the estimated effects on the H3K27ac levels.
pdat5 <- data.frame(effect = obj_plot$effect_s[1,],
                    up     = obj_plot$effect_s[2,],
                    low    = obj_plot$effect_s[3,],
                    pos    = obj_plot$pos_H3Kac_effect)
n     <- length(obj_plot$peak_pos)
pdat5 <- interpolate_effect_estimates(pdat5,obj_plot$peak_pos[seq(2,n-1)])
pdat5 <- subset(pdat5,
                pos >= pos0 &
                pos <= pos1)
pdat5 <- transform(pdat5,
                   pos    = pos/1e6,
                   effect = -effect,
                   low    = -up,
                   up     = -low,
                   y      = 0)
pdat5 <- list(zero_effects = subset(pdat5,low <= 0 & up >= 0),
              nonzero_effects = subset(pdat5,!(low <= 0 & up >= 0)))
p5 <- ggplot() +
  geom_hline(yintercept = 0,linetype = "dotted") +
  geom_linerange(data = pdat5$nonzero_effects,
                 mapping = aes(x = pos,ymin = low,ymax = up),
                 color = "dodgerblue") +
  geom_point(data = pdat5$nonzero_effects,
             mapping = aes(x = pos,y = effect),,
             color = "dodgerblue",size = 0.75) +
  geom_point(data = pdat5$zero_effects,
             mapping = aes(x = pos,y = y),
             color = "dodgerblue",size = 0.75) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(207,208,0.05),
                     labels = NULL) +
  labs(x = "",y = "effect") +
  theme_cowplot(font_size = 9)

# The sixth panel shows the raw data.
pdat6 <- obj_plot$count_df
names(pdat6) <- c("AA","AG","GG","pos")
pdat6 <- subset(pdat6,pos >= pos0 & pos <= pos1)
rows1 <- with(pdat6,which(pmax(AA,AG,GG) >= 5))
rows2 <- with(pdat6,which(pmax(AA,AG,GG) < 5))
rows2 <- sample(rows2,5000)
rows  <- c(rows1,rows2)
pdat6 <- pdat6[rows,]
pdat6 <- transform(pdat6,pos = pos/1e6)
pdat6 <- melt(pdat6,id.vars = "pos",variable.name = "genotype",
              value.name = "count")
pdat6 <- transform(pdat6,genotype = factor(genotype))
rows  <- order(pdat6$genotype,decreasing = TRUE)
pdat6 <- pdat6[rows,]
p6 <- ggplot(pdat6,aes(x = pos,y = count,color = genotype)) +
  geom_point(size = 0.35) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  scale_color_manual(values = c("darkblue","darkviolet","darkorange")) +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(207,208,0.05),
                     labels = NULL) +
  ylim(2,24) +
  labs(x = "") + 
  theme_cowplot(font_size = 9)
             
# The seventh panel shows the genes.
pdat7 <- subset(genes,
                chromosome == "chr1" &
                end > pos0 &
                start < pos1)
pdat7 <- transform(pdat7,tss = ifelse(strand == "+",start,end))
pdat7 <- transform(pdat7,
                   start = start/1e6,
                   end   = end/1e6,
                   tss   = tss/1e6)
n <- nrow(pdat7)
pdat7$y <- seq(0,1,length.out = n)
p7 <- ggplot(pdat7,aes(x = start,xend = end,y = y,yend = y,
                       label = gene_name)) +
  geom_segment(color = "dodgerblue",linewidth = 0.5) +
  geom_point(mapping = aes(x = tss),color = "dodgerblue",size = 1.5,
             shape = 18) +
  geom_text(color = "black",size = 2.25,fontface = "italic",
            hjust = "right",nudge_x = -0.003) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(207,208,0.05)) +
  scale_y_continuous(limits = c(-0.1,1.1),breaks = NULL) +
  labs(x = "base-pair position on chromosome 1 (Mb)",y = "") + 
  theme_cowplot(font_size = 9)

# Save the full figure to a PDF.
print(plot_grid(p1,p2,p3,p4,p5,p6,p7,nrow = 7,ncol = 1,align = "v"))
ggsave("CR1CR2_zoomin_plot.pdf",
       plot_grid(p1,p2,p3,p4,p5,p6,p7,nrow = 7,ncol = 1,align = "v"),
       height = 5.75,width = 4.5)

stop()

data_track =plot_df[ which(  plot_df$study ==study & plot_df$region==gene_name ),]
 
data_track_CS1 =plot_df [ which(  plot_df$study == study & plot_df$CS1& plot_df$region==gene_name  ),]  

# The third panel shows the CR2 eQTL p-values.

library(ggplot2)
library(tidyr)
library(AnnotationHub)
library(org.Hs.eg.db)
library(GenomicRanges)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(fsusieR)
library(dplyr) 

cex=0.6
chr=1
study="AD_Bellenguez_2022"
 #### AD   -----

plot_df= obj_plot$ plot_df 
view_win= obj_plot$ view_win

data_track = plot_df [ which(   plot_df $study == study   ),]
 
data_track_CS1 =plot_df [ which(  plot_df$study == study & plot_df$CS1   ),]   


t1= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track$pos , end = data_track$pos )),
                data = matrix(data_track$`-log10(P)` , nrow=1), genome = "hg19",
                ylim =c( min(data_track$`-log10(P)`), max(data_track$`-log10(P)`)+0.2),
                type = "p", col = "black",  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
                background.title = "white",name="AD") ) # Change title color to black




t2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track_CS1$pos   , end = data_track_CS1$pos  )),
                data = matrix(data_track_CS1$`-log10(P)`  , nrow=1), genome = "hg19", 
                ylim =c( min(data_track$`-log10(P)`), max(data_track$`-log10(P)`)+0.2),
                type = "p", col = "black", cex=1.5,
                fill=  "royalblue",
                pch=c(25 ),# Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
                background.title = "white",name="AD") ) # Change title color to black


t3= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track_CS1$pos [c(2,4,5)]   ,
                                                                 end = data_track_CS1$pos [c(2,4,5)]  )),
                data = matrix(data_track_CS1$`-log10(P)`[c(2,4,5)]  , nrow=1), genome = "hg19", 
                ylim =c( min(data_track$`-log10(P)`), max(data_track$`-log10(P)`)+0.2),
                type = "p", col = "red", cex=1.5,
                fill=  "royalblue",
                pch=c(25 ),# Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
                background.title = "white",name="AD") ) # Change title color to black

otAD <- OverlayTrack(trackList=list(    t1,  t2,t3  ),
                     background.title = "white")


AD_cs =data_track_CS1
plotTracks( otAD )




#### CR1 panel ---- 


gene_name="CR1"
study="DLPFC_DeJager_eQTL" 
data_track =plot_df[ which(  plot_df$study ==study & plot_df$region==gene_name ),]
 
data_track_CS1 =plot_df [ which(  plot_df$study == study & plot_df$CS1& plot_df$region==gene_name  ),]  


t1= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track$pos , end = data_track$pos )),
                data = matrix(data_track$`-log10(P)` , nrow=1), genome = "hg38", 
                type = "p", col = "black",  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
                background.title = "white",name="CR1") ) # Change title color to black




t2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track_CS1$pos , end = data_track_CS1$pos )),
                data = matrix(data_track_CS1$`-log10(P)` , nrow=1), genome = "hg38", 
                ylim =c( min(data_track$`-log10(P)`), max(data_track$`-log10(P)`)),
                type = "p",  col = "black", cex=1.5,
                fill=  "royalblue",
                pch=c(25 ),
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
                background.title = "white",name="CR1") ) # Change title color to black


otCR1 <- OverlayTrack(trackList=list(    t1, t2 ),
                      background.title = "white")

CR1_cs=data_track_CS1
plotTracks( otCR1 )




#### CR2 panel ---- 


gene_name="CR2"
study="DLPFC_DeJager_eQTL" 
data_track =plot_df[ which(  plot_df$study ==study & plot_df$region==gene_name ),]
 
data_track_CS1 =plot_df [ which(  plot_df$study == study & plot_df$CS1& plot_df$region==gene_name  ),]  #%>% filter(study == "DLPFC_DeJager_eQTL", CS1)#"maroon"



t1= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track$pos , end = data_track$pos )),
                data = matrix(data_track$`-log10(P)` , nrow=1), genome = "hg38", 
                type = "p", col = "black",  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
                background.title = "white",name="CR2") ) # Change title color to black




t2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track_CS1$pos , end = data_track_CS1$pos )),
                data = matrix(data_track_CS1$`-log10(P)` , nrow=1), genome = "hg38", 
                ylim =c( min(data_track$`-log10(P)`), max(data_track$`-log10(P)`)),
                type = "p",  col = "black", cex=1.5,
                fill=  "royalblue",
                pch=c(25 ),
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
                background.title = "white",name="CR2") ) # Change title color to black


otCR2 <- OverlayTrack(trackList=list(    t1, t2 ),
                      background.title = "white")

CR2_cs=data_track_CS1
plotTracks( otCR2 )




##### effect plot  -----

effect_s=obj_plot$effect_s
 
positions=obj_plot$pos_H3Kac_effect
pos= obj_plot$peak_pos
chrom=1
plot_list=list()
df_list=list()
widthtick=2500


#here I need to redo the error bar my self
#I create all the plot separately for each error bar
#then I stack then for the one I am interested in
for ( i in 2:(length(pos)-1))
  
{
  
  
  up =  which(positions>=  pos[i] )
  low = which( pos[i] >=positions  ) 
  
  idx=   ifelse(  abs( pos[i]-positions[min(up)])<  abs( pos[i]-positions[ max(low)]),
                  min(up),
                  max(low))[1]
  
  
  du = abs( pos[i]-positions[min(up)])
  di = abs( pos[i]-positions[ max(low)])
  dupdi=du+di
  
  effect =  (1 - du/dupdi ) *effect_s[1, min(up)]+ (1 - di/dupdi )* effect_s[1, max(low)]
  ci_upper =   (1 - du/dupdi ) *effect_s[2, min(up)]+ (1 - di/dupdi ) *effect_s[2, max(low)]  
  ci_lower =  (1 - du/dupdi )* effect_s[3, min(up)]+ (1 - di/dupdi )* effect_s[3, max(low)] 
  
  if(ci_lower>0 | ci_upper<0){
    df_list[[i-1]]= data.frame(effect=effect,
                               ci_lower=ci_lower,
                               ci_upper=ci_upper,
                               pos=pos[i])
    
    # DataTrack for effect size points (blue dots)
    dTrack_points <- DataTrack(start = pos[i ], end = pos[i ], genome = "hg38", chromosome = chrom,
                               data = effect, name = "Effect Size", type = "p",
                               ylim =c( min( c(effect_s)),max(c(effect_s)  )),
                               col = "royalblue", pch = 16, cex =   0.41,
                               background.title = "white" )
    
    # Create GRanges for vertical error bars
    gr_errors <- GRanges(seqnames = chrom,
                         ranges = IRanges(start = pos[i ], end = pos[i ]),
                         lower = ci_lower,
                         upper = ci_upper)
    # AnnotationTrack for vertical error bars
    error_bar_track <- DataTrack(start = pos[i ], end = pos[i ], genome = "hg38", chromosome = chrom,
                                 data = matrix(
                                   c(ci_lower, ci_upper), ncol=1
                                 ), name = "Effect Size", type = "l",
                                 ylim =c( min( c(effect_s)),max(c(effect_s)  )),
                                 col = "royalblue", pch = 16, cex = 1.2,
                                 lwd=2,
                                 ylim =c( min( c(effect_s)),max(c(effect_s)  )),
                                 
                                 background.title = "white"
    )
    
    tick_up <- DataTrack(start = ( pos[i] +-widthtick:widthtick),
                         end = pos[i] +-widthtick:widthtick+1,
                         genome = "hg38", chromosome = chrom,
                         data = matrix(
                           rep( ci_upper, (2*widthtick+1)), 
                           nrow= 1)  ,
                         name = "Effect Size", type = "l",
                         lwd=2,
                         ylim =c( min( c(effect_s)),max(c(effect_s)  )),
                         col = "royalblue", pch = 16, cex = 1.2,
                         ylim =c( min( c(effect_s)),max(c(effect_s)  )),
                         background.title = "white"
    ) 
    
    
    
    tick_low <- DataTrack(start = ( pos[i] +-widthtick:widthtick),
                          end = pos[i] +-widthtick:widthtick+1,
                          genome = "hg38", chromosome = chrom,
                          data = matrix(
                            rep( ci_lower, (2*widthtick+1)), 
                            nrow= 1)  ,
                          lwd=2,
                          name = "Effect Size", type = "l",
                          ylim =c( min( c(effect_s)),max(c(effect_s)  )),
                          col = "royalblue", pch = 16, cex = 1.2,
                          ylim =c( min( c(effect_s)),max(c(effect_s)  ))
    ) 
    
    
    tt = OverlayTrack(trackList = list( tick_low,tick_up,error_bar_track,  dTrack_points))
    #plotTracks(tt)
    
    plot_list[[i-1]] <- OverlayTrack(trackList = tt,
                                     background.title = "white")
  }else{
    
    df_list[[i-1]]= data.frame(effect=effect,
                               ci_lower=ci_lower,
                               ci_upper=ci_upper,
                               pos=pos[i])
    
    
    
    
    plot_list[[i-1]] <-   DataTrack(start = pos[i ], end = pos[i ], genome = "hg38", chromosome = chrom,
                                    data = 0*effect, name = "Effect Size", type = "p",
                                    ylim =c( min( c(effect_s)),max(c(effect_s)  )),
                                    col = "black", pch = 16, cex =  2,
                                    background.title = "white" )
  }
  
  
}



tt= do.call( rbind , df_list)

view_win <- c(207317782, 207895513)
#Keeping the plot within the region of interest
idl= which( pos > view_win[1] & pos < view_win[2])-1
#stack then into a single plot
total_overlay= OverlayTrack( trackList =plot_list[idl],
                             background.title = "white")





effect0=       rep(0,ncol( effect_s ))
group_cred= c( 0)
group_colors <- c("black"  )

group_lwd= c(1)
haQTL_pos0 =   DataTrack(range = GRanges(seqnames = chr,
                                         ranges = IRanges(start = positions,
                                                          end = positions + 1)),
                         data = effect0, genome = "hg38",
                         groups= group_cred,
                         ylim =c( min( c(effect_s)),max(c(effect_s)  )) ,
                         lwd = group_lwd,
                         rotation.title = 90,
                         name ="effect H3k9ac",
                         type = c(  "l" ),
                         col = group_colors,
                         
                         track.margin = 0.05,
                         cex=1.5,# Use color column from df_plot
                         track.margin = 0.05, # Reduce margin between track and title
                         cex.title = 0.6,     # Reduce title size
                         cex.axis = 0.6,      # Reduce axis text size
                         col.axis = "black",  # Change axis color to black
                         col.title = "black",
                         background.title = "white",
                         legend = FALSE  # Remove legend
)

 #add a line a 0 


fsusie_ha_plot <- OverlayTrack(trackList = list( haQTL_pos0,
                                                 OverlayTrack(trackList =plot_list[idl])
),
background.title = "white"
)

plotTracks(fsusie_ha_plot  )



##### PIP plots  ----- 
pip_df= obj_plot$pip_df 



col_names =obj_plot$name_SNP 
pos_SNP_HA <-  as.numeric(gsub("chr[0-9XY]+\\.([0-9]+)\\..*", "\\1", col_names))
pip_df=obj_plot$pip_df


data_ha =pip_df[which( pip_df$study =="ROSMAP_DLPFC_haQTL"&pip_df$cs_coverage_0.95_min_corr==5  ),]
tdf     = plot_df[ which(  plot_df$study ==study & plot_df$region==gene_name ),]
tdf$z   =tdf$z*0 

 
# subset and shfiting data so it look nice
t_dat=obj_plot$pip_fsusie_obj
t_dat[obj_plot$cs_fsusie_obj[[5]]]=t_dat[obj_plot$cs_fsusie_obj[[5]]]+0.05

t_0= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = pos_SNP_HA  , end = pos_SNP_HA  )),
                 data = matrix(t_dat, nrow=1), genome = "hg38", 
                 ylim =c( 0, 0.5),
                 type = "p", col = "black",
                 # cex=1.5,# Use color column from df_plot
                 track.margin = 0.05, # Reduce margin between track and title
                 cex.title = 0.6,     # Reduce title size
                 cex.axis = 0.6,      # Reduce axis text size
                 col.axis = "black",  # Change axis color to black
                 col.title = "black",rotation.title = 90,cex.title = cex,
                 background.title = "white",name="PIP \n H3K9ac") )  
t_ha0= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_ha$pos , end = data_ha$pos )),
                   data = matrix(data_ha$pip+0.05 , nrow=1), genome = "hg38", 
                   #offset for the plot to visualize
                   ylim =c( 0, 0.5),
                   type = "p",  
                   
                   col = "black", cex=1.5,
                   fill=  "royalblue",
                   pch=c(25 ),
                   
                   cex=1.5,# Use color column from df_plot
                   track.margin = 0.05, # Reduce margin between track and title
                   cex.title = 0.6,     # Reduce title size
                   cex.axis = 0.6,      # Reduce axis text size
                   col.axis = "black",  # Change axis color to black
                   col.title = "black",rotation.title = 90,cex.title = cex,
                   background.title = "white",name="PIP \n H3K9ac") ) # Change title color to black


HA_cs=data_ha

tidx= which( HA_cs$pos %in%AD_cs$pos )
 
t_ha1= ( DataTrack(range = GRanges(seqnames = chr, 
                                   ranges = IRanges(start = data_ha$pos[tidx] , end = data_ha$pos[ tidx] )),
                   data = matrix(data_ha$pip[tidx]+0.03 , nrow=1), genome = "hg38", 
                   #offset for the plot to visualize
                   ylim =c( 0, 0.5),
                   type = "p",  
                   
                   col = "red", cex=1.5,
                   fill=  "royalblue",
                   pch=c(25 ),
                   
                   cex=1.5,# Use color column from df_plot
                   track.margin = 0.05, # Reduce margin between track and title
                   cex.title = 0.6,     # Reduce title size
                   cex.axis = 0.6,      # Reduce axis text size
                   col.axis = "black",  # Change axis color to black
                   col.title = "black",rotation.title = 90,cex.title = cex,
                   background.title = "white",name="PIP \n H3K9ac") ) # Change title color to black


 



t_ha=OverlayTrack(trackList=list( t_0, t_ha0, t_ha1 ),
                  background.title = "white")
plotTracks( t_ha)


list_track=  list( otAD,
                   otCR1,
                   otCR2,
                   t_ha
)

view_win <- c(207317782, 207895513)


### Count ----- 
 

df =obj_plot$count_df 

# summarize count averaging over bin size for different genotype
bin_size=2000
# Create a bin column
df$bin <- floor(df$obs_pos / bin_size)

# Load dplyr and aggregate
library(dplyr)

binned_df <- df %>%
  group_by(bin) %>%
  summarise(
    mean_func0 = sum(mean_func0, na.rm = TRUE),
    mean_func1 = sum(mean_func1, na.rm = TRUE),
    mean_func2 = sum(mean_func2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(bin_start_pos = bin * bin_size)

# View result
head(binned_df)



plot(binned_df$mean_func2, type="l", lwd=2)
lines(binned_df$mean_func1, col="green", lwd=2)
lines(binned_df$mean_func0, col="red", lwd=2)

binned_df$mean_func2= (binned_df$mean_func2)
binned_df$mean_func1= (binned_df$mean_func1)

binned_df$mean_func0= (binned_df$mean_func0)



library(Gviz)
 #preparing the tack
gr_func0 <- GRanges(
  seqnames = "chr1",  # use real chromosome if you have it
  ranges = IRanges(start = binned_df$bin_start_pos,
                   width = bin_size),
  score = binned_df$mean_func0
)

gr_func1 <- gr_func0
mcols(gr_func1)$score <- binned_df$mean_func1

gr_func2 <- gr_func0
mcols(gr_func2)$score <- binned_df$mean_func2
 


# Create DataTracks
track0 <- DataTrack(gr_func0, 
                    type = "hist",  
                    fill = "turquoise",
                    col.histogram = NA, 
                    ylim = c(0, 25 ),#max(c(binned_df$mean_func0, binned_df$mean_func1, binned_df$mean_func2))),
                    cex=1.5, # Use color column from df_plot
                    track.margin = 0.05, # Reduce margin between track and title
                    cex.title = 0.6,     # Reduce title size
                    cex.axis = 0.6,      # Reduce axis text size
                    col.axis = "black",  # Change axis color to black
                    col.title = "black",rotation.title = 90,cex.title = cex,
                    background.title = "white",name="Observed count")
track1 <- DataTrack(gr_func1, 
                    type = "hist",  
                    col.histogram = NA, 
                    fill = "turquoise",
                    ylim = c(0, 25 ),#max(c(binned_df$mean_func0, binned_df$mean_func1, binned_df$mean_func2))),
                    cex=1.5, # Use color column from df_plot
                    track.margin = 0.05, # Reduce margin between track and title
                    cex.title = 0.6,     # Reduce title size
                    cex.axis = 0.6,      # Reduce axis text size
                    col.axis = "black",  # Change axis color to black
                    col.title = "black",rotation.title = 90,cex.title = cex,
                    background.title = "white",name="Observed count")
track2 <- DataTrack(gr_func2, 
                    type = "hist", 
                    col.histogram = NA, 
                    fill = "royalblue" ,
                    ylim = c(0, 25 ),#max(c(binned_df$mean_func0, binned_df$mean_func1, binned_df$mean_func2))),
                    cex=1.5, # Use color column from df_plot
                    track.margin = 0.05, # Reduce margin between track and title
                    cex.title = 0.6,     # Reduce title size
                    cex.axis = 0.6,      # Reduce axis text size
                    col.axis = "black",  # Change axis color to black
                    col.title = "black",rotation.title = 90,cex.title = cex,
                    background.title = "white",name="Observed count")

# Plot
ot_count =OverlayTrack(trackList =  list(  track2,track1 ),track.margin = 0.05,
                       background.title = "white")
 

list_track=  list( otAD,
                   otCR1,
                   otCR2,
                   t_ha,
                   ot_count
)

plotTracks(list_track,
           from = min( plot_df$pos[which(plot_df$study=="AD_Bellenguez_2022")]),
           to=max( plot_df$pos[which(plot_df$study=="AD_Bellenguez_2022")]) )

plotTracks(list_track,
           from =view_win[1],
           to=view_win[2])




### gene track -----


library(biomaRt)
library(GenomicRanges)

# Define the genomic region
chr <- "chr1"
start_pos <-  min( plot_df$pos[which( 
  plot_df$study=="DLPFC_DeJager_eQTL")])
end_pos <-max( plot_df$pos[which(plot_df$study=="AD_Bellenguez_2022")])

# Use biomaRt to fetch gene annotations from Ensembl
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# Get gene and transcript information
genes <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", 
                 "strand", "ensembl_gene_id", "ensembl_transcript_id", "external_gene_name"),
  filters = c("chromosome_name", "start", "end"),
  values = list("1", start_pos, end_pos),
  mart = mart
)

# Get exon-level information
exons <- getBM(
  attributes = c("chromosome_name", "exon_chrom_start", "exon_chrom_end",
                 "strand", "ensembl_gene_id", "ensembl_transcript_id", 
                 "ensembl_exon_id", "external_gene_name"),
  filters = c("chromosome_name", "start", "end"),
  values = list("1", start_pos, end_pos),
  mart = mart
)

# Check if any genes were retrieved
if (nrow(genes) == 0) {
  stop("No gene data retrieved. Check chromosome and coordinates.")
}

# Ensure strand is correctly formatted
exons$strand <- ifelse(exons$strand == 1, "+", "-")

# Keep only one isoform per gene (longest transcript)
genes <- genes[order(genes$external_gene_name, genes$end_position - genes$start_position, decreasing = TRUE), ]
genes <- genes[!duplicated(genes$external_gene_name), ]

# Filter exons to match selected transcripts
exons <- exons[exons$ensembl_transcript_id %in% genes$ensembl_transcript_id, ]

# Rename columns to match GeneRegionTrack expectations
exons <- exons[, c("chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand", "ensembl_gene_id", "ensembl_transcript_id", "ensembl_exon_id", "external_gene_name")]
colnames(exons) <- c("chromosome", "start", "end", "strand", "gene", "transcript", "exon", "symbol")

# Convert to GeneRegionTrack-compatible data frame
exons_df <- data.frame(
  chromosome = paste0("chr", exons$chromosome),
  start = exons$start,
  end = exons$end,
  strand = exons$strand,
  gene = exons$gene,
  transcript = exons$transcript,
  exon = exons$exon,
  symbol = exons$symbol,
  feature = "exon"  # Mark exons so GViz can differentiate introns automatically
)

# Create the GeneRegionTrack with exon/intron display
gene_track <- GeneRegionTrack(
  exons_df,
  genome = "hg38",
  chromosome = chr,
  start = start_pos,
  end = end_pos,
  name = "Genes",
  showId = TRUE,
  transcriptAnnotation = "symbol",
  col = "black",
  fill = "blue",
  
  col.axis = "black",col.title = "black",
  rotation.title = 90,cex.title = cex,
  col = "salmon",fill = "salmon",
  background.title = "white"
)
gtrack <- GenomeAxisTrack()




list_track=  list(     
  otAD,
  otCR1,
  otCR2  ,
  
  t_ha,
  fsusie_ha_plot ,
  ot_count,
  gene_track,
  gtrack 
  
)

 
plotTracks(list_track,
           from =view_win[1],
           to=view_win[2])

 



file_path <- file.path(folder_path, "CR1_CR2.pdf")
pdf(file_path, width = 8.27, height = 11.69)  # A4 in inches


plotTracks(list_track,
           from =207441000 ,   #view_win[1],
           to=   207683000  ,  #view_win[2]  ,
           frame = TRUE ,
           sizes = c(0.75,0.75, 0.75,0.35,  0.5,0.75, 0.5,0.3),
           #fontsize  = 15
           cex.main=1.2, cex.title = 1.
)

grid.text(
  "rs679515",
  x = 0.65,
  y =0.46,
  gp = gpar(col = "black", fontsize = 10)
)
grid.text(
  "rs10863418",
  x = 0.75,
  y =0.43 ,
  gp = gpar(col = "black", fontsize = 10)
)
grid.text(
  "rs4844610",
  x = 0.84,
  y =0.49,
  gp = gpar(col = "black", fontsize = 10)
)

dev.off()



####  global PIP plot ---- 




plot_colors <-          c("black" , "steelblue4", 
                          "green4", "deeppink1",
                          "#6A3D9A","royalblue",  
                          "darkturquoise", "green1",
                          "yellow4")
SNP_in_cs = c(unlist( obj_plot$cs_fsusie_obj )) 

idx= which(obj_plot$pip_fsusie_obj  >0.05  )
obj_plot$pip_fsusie_obj [idx[ -which(idx %in% SNP_in_cs  )] ]=0.0
pos_SNP =  pos_SNP_HA
 
point_size = 1.25
L <- length(obj_plot$cs_fsusie_obj)
y <- obj_plot$pip_fsusie_obj
font_size = 10
col_y <- rep(0,length(y))
for (l in 1:L) {
  col_y[which(1:length(y) %in% obj_plot$cs_fsusie_obj[[l]])] <- l
}
point_shape <- rep(19,length(pos_SNP))
CS <- factor(col_y,
             levels = 0:L,
             labels = c("none",1:L))
df <- data.frame(y = y,CS = CS,pos_SNP)
P_pip_ha = ggplot(df,aes(y = y,x = pos_SNP,color = CS)) +
  geom_point(size =2.25,
             shape = point_shape) +
  scale_color_manual("CS",values = plot_colors) +
  labs(x = "SNP",y = "PIP" )  +
  theme_minimal( )+
  theme(legend.position="none")+
  
  geom_rect(aes(xmin = view_win[1], 
                xmax = view_win[2],
                ymin = -0.01,
                ymax = 1.02), 
            alpha = 0.0, color = "red")



file_path <- file.path(folder_path, "PIP_ha.pdf")
pdf(file_path, width =  11.69, height =8.27)  # A4 in inches


P_pip_ha  


dev.off() 
