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
 
pdat =fread(paste0(path,"/data/fig_4_data/GWAS_eQTL_sumstat.chr6_44880827_48309905.tsv"))

dim(pdat)
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





res$







#### meqtl -----



res <- readRDS("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/CD2AP/ROSMAP_mQTL.chr6_44880827_48309905.fsusie_mixture_normal_top_pc_weights.input_data.rds")

fsusie_obj_me = res$`chr6:44880827-48309905`$ROSMAP_DLPFC_mQTL$fsusie_result
rm(res)

res_me <- readRDS("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/ROSMAP_mQTL.chr6_44880827_48309905.fsusie_mixture_normal_top_pc_weights.input_data.rds")
## to work from here
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

## reprendre de la WW---- 
plotTracks(list_track)

plot

Y= as.data.frame(res_me$residual_Y)

col_names <-        colnames(as.data.frame(res_me$X_data)) 

pos_SNP_me <-  as.numeric(gsub("chr[0-9XY]+\\.([0-9]+)\\..*", "\\1", col_names))


X=as.data.frame(res_me$residual_X)
pos = as.data.frame(res_me$Y_coordinates) #use start
pos= pos$start


map_data <- fsusieR:::remap_data(Y=Y,
                                 pos=pos,
                                 
                                 max_scale=10)

outing_grid <- map_data$outing_grid
Y_w= map_data$Y
