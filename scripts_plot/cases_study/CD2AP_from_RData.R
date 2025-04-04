# Create a "zoom-in" plot for the CD2AP locus.
#
# NOTE: download CD2AP_obj.RData and ROSMAP_mQTL.chr6_44880827_48309905.
# fsusie_mixture_normal_top_pc_weights.input_data_1.rds from
# https://uchicago.box.com/s/tt1vgg7vqayfthg0vsbl8gw0sdi0f1uo
library(data.table)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(cowplot)
source("get_gene_annotations.R")
source("interpolate_effect_estimates.R")
gene_file <-
  file.path("../../data/genome_annotations",
    "Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz")
genes <- get_gene_annotations(gene_file)
raw_dat <- readRDS("../../outputs/ROSMAP_mQTL.chr6_44880827_48309905.fsusie_mixture_normal_top_pc_weights.input_data_1.rds")
load("../../outputs/CD2AP_obj.RData")

pos0 <- 47.375e6
pos1 <- 47.725e6
key_marker <- 47472829/1e6

# The top panel shows the Alzheimer's Disease (AD) association
# p-values.
pdat1 <- as.data.frame(obj_plot$pdat)
pdat1 <- subset(pdat1,
                study == "AD_Bellenguez_2022" &
                pos >= pos0 &
                pos <= pos1)
pdat1 <- pdat1[c("variant_alternate_id","pos","CS1","-log10(P)")]
pdat1 <- transform(pdat1,pos = pos/1e6)
# > nrow(subset(pdat1,CS))
# 34
names(pdat1) <- c("id","pos","CS","pval")
ids <- pdat1$id
ids[] <- NA
ids[pdat1$id == "chr6:47472829:C:A"] <- "rs9369695"
ids[pdat1$id == "chr6:47615928:C:T"] <- "rs12195738"
pdat1$id <- ids
p1 <- ggplot(pdat1,aes(x = pos,y = pval,color = CS,label = id)) +
  geom_point(size = 0.75) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  geom_text_repel(size = 2.25,color = "dimgray",segment.color = "dimgray",
                  min.segment.length = 0,max.overlaps = Inf) +
  scale_color_manual(values = c("black","dodgerblue")) +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(45,48,0.05),
                     labels = NULL) +
  ylim(0,15) +
  labs(x = "",y = "AD") + 
  theme_cowplot(font_size = 9)

# The second panel shows the CD2Ap eQTL p-values.
pdat2 <- as.data.frame(obj_plot$pdat)
pdat2 <- subset(pdat2,
                region == "CD2AP" &
                study == "Exc_DeJager_eQTL" &
                pos >= pos0 &
                pos <= pos1)
pdat2$CS <- as.numeric(NA)
pdat2[which(pdat2$CS1),"CS"] <- 1
pdat2 <- pdat2[c("variant_alternate_id","pos","CS","-log10(P)")]
names(pdat2) <- c("id","pos","CS","pval")
pdat2 <- transform(pdat2,
                   pos = pos/1e6,
                   CS  = factor(CS))
# > nrow(subset(pdat2,CS == 1))
# 80
ids <- pdat2$id
ids[] <- NA
ids[pdat2$id == "chr6:47472829:C:A"] <- "rs9369695"
pdat2$id <- ids
p2 <- ggplot(pdat2,aes(x = pos,y = pval,color = CS,label = id)) +
  geom_point(size = 0.75) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  scale_color_manual(values = c("limegreen","gold"),na.value = "black") +
  geom_text_repel(size = 2.25,color = "dimgray",segment.color = "dimgray",
                  min.segment.length = 0,max.overlaps = Inf) +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(45,48,0.05),
                     labels = NULL) +
  ylim(0,25) +
  labs(x = "",y = "CD2AP eQTL") + 
  theme_cowplot(font_size = 9)

# The third panel shows the mSNP PIPs.
ids   <- names(obj_plot$fsusie_obj_me$pip)
pdat3 <- data.frame(id  = as.character(NA),
                   pos = sapply(strsplit(ids,":",fixed = TRUE),"[[",2),
                   pip = obj_plot$fsusie_obj_me$pip,
                   cs  = as.character(NA),
                   stringsAsFactors = FALSE)
pdat3 <- transform(pdat3,pos = as.numeric(pos))
n <- length(obj_plot$fsusie_obj_me$sets$cs)
for (i in 1:n) {
  snps <- names(obj_plot$fsusie_obj_me$sets$cs[[i]])
  pdat3[snps,"cs"] <- i
}
pdat3 <- subset(pdat3,pos >= pos0 & pos <= pos1)
pdat3 <- transform(pdat3,pos = pos/1e6)
# > subset(pdat3,!is.na(cs))
#                     id   pos    pip cs
# chr6:47405566:G:A <NA> 47.41 0.5000 16
# chr6:47405711:A:G <NA> 47.41 0.5000 16
# chr6:47472829:C:A <NA> 47.47 0.9645 13
# chr6:47664375:T:C <NA> 47.66 0.3333  8
# chr6:47664664:A:T <NA> 47.66 0.3333  8
# chr6:47665074:C:T <NA> 47.67 0.3333  8
pdat3[c("chr6:47405566:G:A",
        "chr6:47472829:C:A",
        "chr6:47664375:T:C"),"id"] <-
    c("rs4715009","rs9369695","rs5013495")
p3 <- ggplot(pdat3,aes(x = pos,y = pip,color = cs,label = id)) +
  geom_point(size = 0.75) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  geom_text_repel(size = 2.25,color = "dimgray",segment.color = "dimgray",
                  min.segment.length = 0,max.overlaps = Inf) +
  scale_color_manual(values = c("magenta","darkorange","darkviolet"),
                     na.value = "black") +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(45,48,0.05),
                     labels = NULL) +
  labs(x = "",y = "mSNP PIP") + 
  theme_cowplot(font_size = 9)

# The fourth panel shows the estimated effects on the methylation levels.
pdat4 <- data.frame(effect = obj_plot$effect_s[1,],
                    up     = obj_plot$effect_s[2,],
                    low    = obj_plot$effect_s[3,],
                    pos    = obj_plot$pos_est_effect)
n     <- length(obj_plot$me_pos)
pdat4 <- interpolate_effect_estimates(pdat4,obj_plot$me_pos[seq(2,n-1)])
pdat4 <- subset(pdat4,
                pos >= pos0 &
                pos <= pos1)
pdat4 <- transform(pdat4,
                   pos    = pos/1e6,
                   effect = -effect,
                   low    = -up,
                   up     = -low,
                   y      = 0)
pdat4 <- list(zero_effects    = subset(pdat4,low <= 0 & up >= 0),
              nonzero_effects = subset(pdat4,!(low <= 0 & up >= 0)))
p4 <- ggplot() +
  geom_hline(yintercept = 0,linetype = "dotted") +
  geom_linerange(data = pdat4$nonzero_effects,
                 mapping = aes(x = pos,ymin = low,ymax = up),
                 color = "dodgerblue") +
  geom_point(data = pdat4$nonzero_effects,
             mapping = aes(x = pos,y = effect),,
             color = "dodgerblue",size = 0.75) +
  geom_point(data = pdat4$zero_effects,
             mapping = aes(x = pos,y = y),
             color = "dodgerblue",size = 0.75) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(45,48,0.05),
                     labels = NULL) +
  labs(x = "",y = "effect of A allele") +
  theme_cowplot(font_size = 9)

# The fifth panel shows the raw data.
## pdat5 <- obj_plot$count_df
## names(pdat5) <- c("AA","AG","GG","pos")
## pdat5 <- subset(pdat6,pos >= pos0 & pos <= pos1)
## rows1 <- with(pdat6,which(pmax(AA,AG,GG) >= 5))
## rows2 <- with(pdat6,which(pmax(AA,AG,GG) < 5))
## rows2 <- sample(rows2,5000)
## rows  <- c(rows1,rows2)
## pdat5 <- pdat6[rows,]
## pdat5 <- transform(pdat6,pos = pos/1e6)
## pdat5 <- melt(pdat6,id.vars = "pos",variable.name = "genotype",
##               value.name = "count")
## pdat5 <- transform(pdat6,genotype = factor(genotype))
## rows  <- order(pdat6$genotype,decreasing = TRUE)
## pdat5 <- pdat6[rows,]
## p6 <- ggplot(pdat6,aes(x = pos,y = count,color = genotype)) +
##   geom_point(size = 0.35) +
##   # geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
##   scale_color_manual(values = c("darkblue","darkviolet","darkorange")) +
##   scale_x_continuous(limits = c(pos0,pos1)/1e6,
##                      breaks = seq(207,208,0.05),
##                      labels = NULL) +
##   ylim(2,24) +
##   labs(x = "") + 
##   theme_cowplot(font_size = 9)

# The sixth panel shows the genes.
pdat6 <- subset(genes,
                chromosome == "chr6" &
                end > pos0 &
                start < pos1)
pdat6 <- transform(pdat6,tss = ifelse(strand == "+",start,end))
pdat6 <- transform(pdat6,
                   start = start/1e6,
                   end   = end/1e6,
                   tss   = tss/1e6)
n <- nrow(pdat6)
pdat6$y <- seq(0,1,length.out = n)
p6 <- ggplot(pdat6,aes(x = start,xend = end,y = y,yend = y,
                       label = gene_name)) +
  geom_segment(color = "dodgerblue",linewidth = 0.5) +
  geom_point(mapping = aes(x = tss),color = "dodgerblue",size = 1.5,
             shape = 18) +
  geom_text(color = "black",size = 2.25,fontface = "italic",
            hjust = "right",nudge_x = -0.003) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(45,48,0.05)) +
  scale_y_continuous(limits = c(-0.1,1.1),breaks = NULL) +
  labs(x = "base-pair position on chromosome 6 (Mb)",y = "") + 
  theme_cowplot(font_size = 9)

# Save the full figure to a PDF.
print(plot_grid(p1,p2,p3,p4,p6,nrow = 5,ncol = 1,align = "v"))
# TO DO.

stop()

library(ggplot2)
library(tidyr)
library(AnnotationHub)
library(org.Hs.eg.db)
library(GenomicRanges)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(fsusieR)
library(dplyr)
library(data.table)
folder_path=  paste0(getwd(),
                     "/plot/"
)
load( paste0(folder_path,"CD2AP_obj.RData"))


dat=obj_plot$ dat
pdat=obj_plot$ pdat 

pdat2=obj_plot$ pdat2

extract_snp_position <- function(snp_string) {
  # Split the input string by ':'
  parts <- unlist(strsplit(snp_string, ":"))
  
  # Extract the position
  position <- as.numeric(parts[2])
  
  return(position)
}
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

cex=0.6

##" cas 1 ---- 

 
view_win <- c(47.3e6, 47.75e6)  

view_win=obj_plot$view_win

chr =  paste("chr",6, sep = "")





### AD GWAS panel -----

study="AD_Bellenguez_2022"
idx=which (pdat$study == study )
pdat1= pdat[idx, ]
dim(pdat1)
#pdat1$X.log10.P. 
t1= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start =pdat1$pos ,
                                                                 end = pdat1$pos )),
                data = matrix(-log10(pdat1$pvalue) , nrow=1), genome = "hg19",
                ylim =c( min(-log10(pdat1$pvalue)), max(-log10(pdat1$pvalue))+0.2),
                type = "p", col = "black",  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
                background.title = "white",name="AD") ) # Change title color to black

plotTracks(t1)

pdat1CS = pdat1[which(pdat1$CS1),]
t2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = pdat1CS$pos   , 
                                                                 end = pdat1CS$pos   )),
                data = matrix(-log10(pdat1CS$pvalue) , nrow=1), genome = "hg19", 
                ylim =c( min(-log10(pdat1$pvalue)), max(-log10(pdat1$pvalue))+0.2),
                type = "p", col = "black", cex=1.5,
                fill=  "royalblue",
                pch=c(24 ),# Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
                background.title = "white",name="AD") ) # Change title color to black
plotTracks(t2)



otAD <- OverlayTrack(trackList=list(    t1,  t2  ),
                     background.title = "white")



plotTracks( otAD , from= view_win[1],
            to= view_win[2])

## CD2AP panel ----- 




pdat2$X.log10.P. 
t1= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start =pdat2$pos ,
                                                                 end = pdat2$pos )),
                data = matrix(pdat2$X.log10.P. , nrow=1), genome = "hg19",
                ylim =c( min(pdat2$X.log10.P.), max(pdat2$X.log10.P.)+0.2),
                type = "p", col = "black",  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
                background.title = "white",name="CD2AP") ) # Change title color to black

plotTracks(t1)

pdat2CS = pdat2[which(pdat2$CS1),]
t2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = pdat2CS$pos   , 
                                                                 end = pdat2CS$pos   )),
                data = matrix(pdat2CS$X.log10.P.  , nrow=1), genome = "hg19", 
                ylim =c( min(pdat2$X.log10.P.), max(pdat2$X.log10.P.)+0.2),
                type = "p", col = "black", cex=1.5,
                fill=  "royalblue",
                pch=c(24 ),# Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
                background.title = "white",name="AD") ) # Change title color to black


oteqTL <- OverlayTrack(trackList=list(    t1,  t2  ),
                       background.title = "white")

plotTracks( oteqTL, from= view_win[1],
            to= view_win[2])








# pip plot  -----








#### meqtl -----
chr=6
cex=1

 
fsusie_obj_me= obj_plot$fsusie_obj_me
 


snp_names=attr( fsusie_obj_me$pip, "names")
pos_SNP_me <- as.numeric(sub("chr[0-9XY]+:([0-9]+):.*", "\\1", snp_names))



#pip_df %>% filter(study %in% c("ROSMAP_DLPFC_mQTL", ""), cs_coverage_0.95 == 7) 
t_me= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = pos_SNP_me , end = pos_SNP_me )),
                  data = matrix(fsusie_obj_me$pip , nrow=1), genome = "hg38", 
                  ylim =c( 0. , 1 ),
                  type = "p", col = "black",
                  
                  cex=1 ,# Use color column from df_plot
                  track.margin = 0.05, # Reduce margin between track and title
                  cex.title = 0.6,     # Reduce title size
                  cex.axis = 0.6,      # Reduce axis text size
                  col.axis = "black",  # Change axis color to black
                  col.title = "black",rotation.title = 90,cex.title = cex,
                  background.title = "white",name="PIP \n DNAm") )


list_cs_plot=list()

list_cs_plot[[1]]=t_me
for ( l in 1 :length(fsusie_obj_me$cs)){
  
  idx=fsusie_obj_me$cs[[l]]
  
  list_cs_plot[[l+1]]= ( DataTrack(range = GRanges(seqnames = chr,
                                                   ranges = IRanges(start = pos_SNP_me[idx] ,
                                                                    end = pos_SNP_me[idx] )),
                                   data = matrix(fsusie_obj_me$pip[idx] , nrow=1), genome = "hg38", 
                                   ylim =c( 0. , 1 ),
                                   type = "p", col = l+1,
                                   
                                   cex=1.5,# Use color column from df_plot
                                   track.margin = 0.05, # Reduce margin between track and title
                                   cex.title = 0.6,     # Reduce title size
                                   cex.axis = 0.6,      # Reduce axis text size
                                   col.axis = "black",  # Change axis color to black
                                   col.title = "black",rotation.title = 90,cex.title = cex,
                                   background.title = "white",name="PIP \n DNAm") )
}



pip_overlay= OverlayTrack( trackList =list_cs_plot,
                           background.title = "white")

plotTracks(pip_overlay)

list_track=  list( otAD,
                   oteqTL,
                   pip_overlay 
)
# we actually care about CS 13
#plot_susiF(fsusie_obj_me)
fsusie_obj_me$cs[[13]]

## reprendre de la WW---- 
plotTracks(list_track)
 

 
positions=obj_plot$pos_est_effect
pos= obj_plot$me_pos
effect_s=obj_plot$effect_s

chrom=6
plot_list=list()
df_list=list()
widthtick=2500
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
                                    col = "black", pch = 16, cex =  .7,
                                    background.title = "white" )
  }
  
  
}



tt= do.call( rbind , df_list)

idl=  1:(length(pos) ) #which( pos > view_win[1] & pos < view_win[2])-1
total_overlay= OverlayTrack( trackList =plot_list ,
                             background.title = "white")










effect0=       rep(0,ncol(effect_s  ))
group_cred= c( 0)
group_colors <- c("black"  )

group_lwd= c(1)
meQTL_pos0 =   DataTrack(range = GRanges(seqnames = chr,
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

fsusie_me_plot <- OverlayTrack(trackList = list(meQTL_pos0,
                                                OverlayTrack(trackList =plot_list )
),
background.title = "white"
)

plotTracks(fsusie_me_plot  )





list_track=  list( otAD,
                   oteqTL,
                   pip_overlay ,
                   fsusie_me_plot 
)
plotTracks(list_track)



