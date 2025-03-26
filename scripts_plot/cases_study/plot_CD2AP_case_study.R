rm(list=ls(
  
))
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
 
  
  path= getwd()
 

dat <- read.table(paste0(  path,"/data/fig_4_data/CD2AP/CD2AP_plot_df.tsv"),header = TRUE,
                  sep = "\t",stringsAsFactors = FALSE)
pdat1 <- subset(dat,
                study == "AD_Bellenguez_2022" &
                  pos >= 47.3e6 & pos <= 47.75e6)
pdat2 <- subset(dat,
                study == "Exc_DeJager_eQTL" &
                  pos >= 47.3e6 & pos <= 47.75e6)

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


#view_win <- c(56300000, 56700000) 
 view_win <- c(47.3e6, 47.75e6)  



chr =  paste("chr",6, sep = "")



47472829
 



### AD GWAS panel -----

 
pdat1$X.log10.P. 
t1= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start =pdat1$pos ,
                                                                 end = pdat1$pos )),
                data = matrix(pdat1$X.log10.P. , nrow=1), genome = "hg19",
                ylim =c( min(pdat1$X.log10.P.), max(pdat1$X.log10.P.)+0.2),
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
                data = matrix(pdat1CS$X.log10.P.  , nrow=1), genome = "hg19", 
                ylim =c( min(pdat1$X.log10.P.), max(pdat1$X.log10.P.)+0.2),
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
 
 

t4= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track_CS1$pos[tidx] ,
                                                                 end = data_track_CS1$pos[tidx] )),
                data = matrix(data_track_CS1$`-log10(P)`[tidx] , nrow=1), genome = "hg19", 
                ylim =c( min(data_track$`-log10(P)`), max(data_track$`-log10(P)`)+0.2),
                type = "p", col = "red", cex=1.5,
                fill=  "royalblue",
                pch=c(24),# Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
                background.title = "white",name="AD") ) 

otAD <- OverlayTrack(trackList=list(    t1,  t2  ),
                     background.title = "white")



plotTracks( otAD , from= view_win[1],
            to= view_win[2])

## EQTL panel ----- 




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
                background.title = "white",name="AD") ) # Change title color to black

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

res <- readRDS("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/CD2AP/ROSMAP_mQTL.chr6_44880827_48309905.fsusie_mixture_normal_top_pc_weights.input_data.rds")
 