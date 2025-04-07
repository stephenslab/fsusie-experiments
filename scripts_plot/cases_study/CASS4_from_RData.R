# Create a "zoom-in" plot for the CASS4 locus.
#
# NOTE: download CASS4_obj.RData from
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
load("../../outputs/CASS4_obj.RData")

pos0 <- 56.35e6
pos1 <- 56.6e6
key_marker <- 56409008/1e6

# The top panel shows the Alzheimer's Disease (AD) association
# p-values.
lpfromz <- function (z)
  (log(2) + pnorm(-abs(z),log.p = TRUE))/log(10)
pdat1 <- as.data.frame(obj_plot$AD_GWAS)
pdat1 <- transform(pdat1,
                   pval = -lpfromz(z))
pdat1 <- subset(pdat1,
                pos >= pos0 &
                pos <= pos1)
pdat1 <- pdat1[c("variant_alternate_id","pos","CS1","pval")]
pdat1 <- transform(pdat1,pos = pos/1e6)
names(pdat1) <- c("id","pos","CS","pval")
# > nrow(subset(pdat1,CS))
# 13
ids <- pdat1$id
ids[] <- NA
ids[pdat1$id == "chr20:56409008:G:C"] <- "rs1884913"
ids[pdat1$id == "chr20:56423488:A:G"] <- "rs6014724"
pdat1$id <- ids
p1 <- ggplot(pdat1,aes(x = pos,y = pval,color = CS,label = id)) +
  geom_point(size = 0.75) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  geom_text_repel(size = 2.25,color = "dimgray",segment.color = "dimgray",
                  min.segment.length = 0,max.overlaps = Inf) +
  scale_color_manual(values = c("black","dodgerblue")) +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(56,57,0.05)) +
  ylim(0,10) +
  labs(x = "",y = "AD") + 
  theme_cowplot(font_size = 9)

# The second panel shows the CASS4 eQTL p-values.
pdat2 <- as.data.frame(obj_plot$qTLdata)
pdat2 <- transform(pdat2,
                   pval = -lpfromz(z))
pdat2 <- subset(pdat2,
                region == "CASS4" &
                study == "DLPFC_DeJager_eQTL" &
                pos >= pos0 &
                pos <= pos1)
pdat2 <- pdat2[c("variant_id","pos","cs_coverage_0.95","pval")]
names(pdat2) <- c("id","pos","CS","pval")
pdat2 <- transform(pdat2,
                   pos = pos/1e6,
                   CS  = factor(CS))
# > nrow(subset(pdat2,CS == 0))
# 11
ids <- pdat2$id
ids[] <- NA
ids[pdat2$id == "chr20:56438160:A:G"] <- "rs6014730"
ids[pdat2$id == "chr20:56409008:G:C"] <- "rs1884913"
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
  # ylim(0,25) +
  labs(x = "",y = "CASS4 eQTL") + 
  theme_cowplot(font_size = 9)

# The third panel shows the mSNP PIPs.
# TO DO.

# Save the full figure to a PDF.
print(plot_grid(p1,p2,nrow = 2,ncol = 1,align = "v"))
# TO DO.

stop()

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
                   pos = pos/1e6,
                   y   = 0)
pdat4 <- list(zero_effects    = subset(pdat4,low <= 0 & up >= 0),
              nonzero_effects = subset(pdat4,!(low <= 0 & up >= 0)))
p4 <- ggplot() +
  geom_hline(yintercept = 0,linetype = "dotted") +
  geom_linerange(data = pdat4$nonzero_effects,
                 mapping = aes(x = pos,ymin = low,ymax = up),
                 color = "dodgerblue") +
  geom_point(data = pdat4$nonzero_effects,
             mapping = aes(x = pos,y = effect),
             color = "dodgerblue",size = 0.75) +
  geom_point(data = pdat4$zero_effects,
             mapping = aes(x = pos,y = y),
             color = "dodgerblue",size = 0.75) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(56,57,0.05),
                     labels = NULL) +
  labs(x = "",y = "effect of C allele") +
  theme_cowplot(font_size = 9)

# The fifth panel shows the genes.
pdat5 <- subset(genes,
                chromosome == "chr6" &
                end > pos0 &
                start < pos1)
pdat5 <- transform(pdat5,tss = ifelse(strand == "+",start,end))
pdat5 <- transform(pdat5,
                   start = start/1e6,
                   end   = end/1e6,
                   tss   = tss/1e6)
n <- nrow(pdat5)
pdat5$y <- seq(0,1,length.out = n)
p5 <- ggplot(pdat5,aes(x = start,xend = end,y = y,yend = y,
                       label = gene_name)) +
  geom_segment(color = "dodgerblue",linewidth = 0.5) +
  geom_point(mapping = aes(x = tss),color = "dodgerblue",size = 1.5,
             shape = 18) +
  geom_text(color = "black",size = 2.25,fontface = "italic",
            hjust = "right",nudge_x = -0.003) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(56,57,0.05)) +
  scale_y_continuous(limits = c(-0.1,1.1),breaks = NULL) +
  labs(x = "base-pair position on chromosome 20 (Mb)",y = "") + 
  theme_cowplot(font_size = 9)

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
path= getwd()
load("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/plot/CASS4_obj.RData")
cex=1
path= getwd()

  AD_GWAS=obj_plot$AD_GWAS
 
qTLdata=obj_plot$qTLdata


chr=20
view_win=c(min(AD_GWAS$pos),max(AD_GWAS$pos))
pdat1= AD_GWAS 
dim(pdat1)
#pdat1$X.log10.P. 

pdat1$lg10p=-pnorm(-abs(pdat1$z), log.p = TRUE) / log(10)

t1= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start =pdat1$pos ,
                                                                 end = pdat1$pos )),
                data = matrix( (pdat1$lg10p) , nrow=1), genome = "hg19",
                ylim =c( min( (pdat1$lg10p)), max( (pdat1$lg10p))+0.2),
                type = "p", col = "black",  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
                background.title = "white",name="AD") ) # Change title color to black

plotTracks(t1)

pdat1CS =pdat1[which ( pdat1$CS1),]  
t2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = pdat1CS$pos   , 
                                                                 end = pdat1CS$pos   )),
                data = matrix( (pdat1CS$lg10p) , nrow=1), genome = "hg19", 
                ylim =c( min( (pdat1$lg10p)), max( (pdat1$lg10p))+0.2),
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









## CASS4 panel ----- 



pdat2 =  qTLdata[which(qTLdata$study=="Mic_DeJager_eQTL"),]
#DLPFC also is consistent
#pdat2 = qTLdata[which(qTLdata$study=="DLPFC_DeJager_eQTL"),]

z_to_neglog10p <- function(z) {
  z <- abs(z)
  log_p <- pnorm(z, lower.tail = FALSE, log.p = TRUE)  # log(1 - Φ(z))
  log10_p <- log_p / log(10)
  -log10(2) - log10_p  # two-sided
} 
pdat2$X.log10.P. <- z_to_neglog10p(pdat2$z)


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
                background.title = "white",name="CASS4") ) # Change title color to black

plotTracks(t1, from= view_win[1],
           to= view_win[2])

pdat2CS = pdat2[which(! is.na( pdat2$cs_coverage_0.95)),]
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

plotTracks(list(otAD, oteqTL),
           
           from = 56407019, to =56433488)




# pip plot  -----










#### meqtl -----



fsusie_obj_me =obj_plot$fsusie_obj_me 

snp_names=attr( fsusie_obj_me$pip, "names")
pos_SNP_me <- as.numeric(sub("chr[0-9XY]+:([0-9]+):.*", "\\1", snp_names))

which(pos_SNP_me %in%   c(56408019,56412160,56413016,56414777,56423488))

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
                  background.title = "white",name="PIP \n H3k4a9ac") )


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
                                   background.title = "white",name="PIP \n H3k4a9ac") )
}



pip_overlay= OverlayTrack( trackList =list_cs_plot,
                           background.title = "white")

plotTracks(pip_overlay)

list_track=  list( otAD,
                   oteqTL,
                   pip_overlay 
)
plotTracks(list_track ,from = 56407019, to =56433488)
# we actually care about CS 13
#plot_susiF(fsusie_obj_me)
fsusie_obj_me$cs[[18]]

 
effect_s= obj_plot$effect_s 
 
positions=obj_plot$pos_est_effect
pos= obj_plot$me_pos
chrom=20
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










effect0=       rep(0,ncol( effect_s ))
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

 


