path= getwd()
data = readRDS(paste0(path , 
                      "/data/fig_4_data/Fig4_data.rds"))

library(ggplot2)
library(tidyr)
library(AnnotationHub)
library(org.Hs.eg.db)
library(GenomicRanges)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(fsusieR)
library(dplyr) 
plot_df <- data$e[[1]]
haQTL_df <- data$e[[2]]
gene_info <- data$e[[3]]
sumstat <- data$e[[4]]
pip_df <- data$e[[5]]


cex=0.6



chr=1
study="AD_Bellenguez_2022"
data_track = plot_df [ which(   plot_df $study == study   ),]
#data_track = data_track[which(    data_track$pos > view_win[1] &  data_track$pos <view_win[2]& plot_df$region==gene_name  ),]

data_track_CS1 =plot_df [ which(  plot_df$study == study & plot_df$CS1   ),]  #%>% filter(study == "DLPFC_DeJager_eQTL", CS1)#"maroon"



t1= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track$pos , end = data_track$pos )),
                data = matrix(data_track$z , nrow=1), genome = "hg19", 
                type = "p", col = "black",  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",
                rotation.title = 0,cex.title = cex,
                background.title = "white",name="AD") ) # Change title color to black




t2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track_CS1$pos , end = data_track_CS1$pos )),
                data = matrix(data_track_CS1$z , nrow=1), genome = "hg19", 
                ylim =c( min(data_track$z), max(data_track$z)),
                type = "p", col = "steelblue", cex=1.5,  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",
                rotation.title = 0,cex.title = cex,
                background.title = "white",name="AD") ) # Change title color to black



otAD <- OverlayTrack(trackList=list(    t1,  t2 ),
                     background.title = "white")



plotTracks( otAD )




#### CR1 panel ---- 


gene_name="CR1"
study="DLPFC_DeJager_eQTL" 
data_track =plot_df[ which(  plot_df$study ==study & plot_df$region==gene_name ),]
#data_track = data_track[which(    data_track$pos > view_win[1] &  data_track$pos <view_win[2]& plot_df$region==gene_name  ),]

data_track_CS1 =plot_df [ which(  plot_df$study == study & plot_df$CS1& plot_df$region==gene_name  ),]  #%>% filter(study == "DLPFC_DeJager_eQTL", CS1)#"maroon"



t1= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track$pos , end = data_track$pos )),
                data = matrix(data_track$z , nrow=1), genome = "hg19", 
                type = "p", col = "black",  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 0,
                background.title = "white",name="CR1") ) # Change title color to black




t2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track_CS1$pos , end = data_track_CS1$pos )),
                data = matrix(data_track_CS1$z , nrow=1), genome = "hg19", 
                ylim =c( min(data_track$z), max(data_track$z)),
                type = "p", col = "steelblue", cex=1.5,  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 0,
                background.title = "white",name="CR1") ) # Change title color to black


otCR1 <- OverlayTrack(trackList=list(    t1, t2 ),
                      background.title = "white")


plotTracks( otCR1 )




#### CR2 panel ---- 


gene_name="CR2"
study="DLPFC_DeJager_eQTL" 
data_track =plot_df[ which(  plot_df$study ==study & plot_df$region==gene_name ),]
#data_track = data_track[which(    data_track$pos > view_win[1] &  data_track$pos <view_win[2]& plot_df$region==gene_name  ),]

data_track_CS1 =plot_df [ which(  plot_df$study == study & plot_df$CS1& plot_df$region==gene_name  ),]  #%>% filter(study == "DLPFC_DeJager_eQTL", CS1)#"maroon"



t1= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track$pos , end = data_track$pos )),
                data = matrix(data_track$z , nrow=1), genome = "hg19", 
                type = "p", col = "black",  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 0,
                background.title = "white",name="CR2") ) # Change title color to black




t2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track_CS1$pos , end = data_track_CS1$pos )),
                data = matrix(data_track_CS1$z , nrow=1), genome = "hg19", 
                ylim =c( min(data_track$z), max(data_track$z)),
                type = "p", col = "steelblue", cex=1.5,  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 0,
                background.title = "white",name="CR2") ) # Change title color to black


otCR2 <- OverlayTrack(trackList=list(    t1, t2 ),
                      background.title = "white")


plotTracks( otCR2 )
##### PIP plots  -----

data_ha =pip_df[which( pip_df$study =="ROSMAP_DLPFC_haQTL"&pip_df$cs_coverage_0.95_min_corr==5  ),]
#  pip_df %>% filter(study %in% c("ROSMAP_DLPFC_haQTL"), cs_coverage_0.95_min_corr == 2)
#data_ha= data_ha[which(data_ha$pos> view_win[1] & data_ha$pos<view_win[2]),]
t_ha= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_ha$pos , end = data_ha$pos )),
                  data = matrix(data_ha$pip , nrow=1), genome = "hg19", 
                  ylim =c( 0, 0.5),
                  type = "p", col = "steelblue",
                  cex=1.5,# Use color column from df_plot
                  track.margin = 0.05, # Reduce margin between track and title
                  cex.title = 0.6,     # Reduce title size
                  cex.axis = 0.6,      # Reduce axis text size
                  col.axis = "black",  # Change axis color to black
                  col.title = "black",rotation.title = 0,cex.title = cex,
                  background.title = "white",name="PIP H3K9ac") ) # Change title color to black


plotTracks( t_ha)


list_track=  list( otAD,
                   otCR1,
                   otCR2,
                   t_ha
)

view_win <- c(207317782, 207895513)
plotTracks(list_track,
           from = min( plot_df$pos[which(plot_df$study=="AD_Bellenguez_2022")]),
           to=max( plot_df$pos[which(plot_df$study=="AD_Bellenguez_2022")]) )

plotTracks(list_track,
           from =view_win[1],
           to=view_win[2])





##### effect plot  -----

res <- readRDS("D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/fsusie_object/ROSMAP_haQTL.chr1_205117782_208795513.fsusie_mixture_normal_top_pc_weights.rds")
fsusie_obj_ha = res$`chr1:205117782-208795513`$ROSMAP_DLPFC_haQTL$fsusie_result
rm(res)



positions = fsusie_obj_ha$outing_grid

idx=5
out=list(effect= fsusie_obj_ha$fitted_func[[idx]],
         cred_band=  fsusie_obj_ha$cred_band[[idx]]  )
effect= out$effect


haQTL_track = DataTrack(range = GRanges(seqnames = chr,
                                        ranges = IRanges(start = positions ,
                                                         end = positions   + 1)),
                        data = effect , genome = "hg38",
                        type = "l", 
                        track.margin = 0.05 ,
                        col.axis = "black",col.title = "black",
                        fontface = "plain",rotation.title = 0,cex.title = cex,
                        background.title = "white",name="effect DNAm")



effect=  out$cred_band[1,]




haQTL_trackcb1  = DataTrack(range = GRanges(seqnames = chr,
                                            ranges = IRanges(start = positions ,
                                                             end = positions + 1)),
                            data = effect , genome = "hg38",
                            type = "l", 
                            track.margin = 0.05 ,lty=2,
                            col.axis = "black",col.title = "black",
                            fontface = "plain",rotation.title = 0,cex.title = cex,
                            background.title = "white",name="effect DNAm")



effect=    out$cred_band[2,]




haQTL_trackcb2  = DataTrack(range = GRanges(seqnames = chr,
                                            ranges = IRanges(start = positions ,
                                                             end = positions + 1)),
                            data = effect , genome = "hg38",
                            type = "l", 
                            track.margin = 0.05 ,lty=2,
                            col.axis = "black",col.title = "black",
                            fontface = "plain",rotation.title = 0,cex.title = cex,
                            background.title = "white",name="effect DNAm")


fsusie_ha_plot <- OverlayTrack(trackList=list( haQTL_track,haQTL_trackcb1,haQTL_trackcb2 ),
                               background.title = "white")
#plotTracks(fsusie_ha_plot  )


list_track=  list( otAD,
                   otCR1,
                   otCR2  ,
                   
                   t_ha,
                   fsusie_ha_plot  
)

view_win <- c(207317782, 207895513)
plotTracks(list_track,
           from = min( plot_df$pos[which(plot_df$study=="AD_Bellenguez_2022")]),
           to=max( plot_df$pos[which(plot_df$study=="AD_Bellenguez_2022")]) )
plotTracks(list_track,
           from =min(positions),
           to=max(positions) )
           
 







## gene track plot ----
library(AnnotationHub)
library(org.Hs.eg.db)
library(GenomicRanges)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Initialize AnnotationHub
ah <- AnnotationHub()

# Load the TxDb for GRCh38
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

chr=paste0("chr",1)
# Extract the relevant genes and exons in the specified region
region_genes <- genes(txdb,columns = c("tx_id","gene_id"))

# Subset the genes and exons to the region of interest.
region_genes <- subsetByOverlaps(region_genes,
                                 GRanges(seqnames = chr,
                                         ranges = IRanges(view_win[1],
                                                          view_win[2])))
cex <- 0.6
# Create a gene region track for the specified region
gene_track <- GeneRegionTrack(txdb,genome = "hg38",chromosome = chr,
                              pos0 = view_win[1],pos1 =view_win[2],name = "",
                              showId = TRUE,geneSymbol = TRUE,
                              col.axis = "black",col.title = "black",
                              transcriptAnnotation = "symbol",
                              rotation.title = 0,cex.title = cex,
                              col = "salmon",fill = "salmon",
                              background.title = "white" )
# Map gene IDs to gene symbols
gene_ids <- unique(unlist(region_genes$gene_id))  # Get unique gene IDs

# Map to gene symbols using org.Hs.eg.db
gene_symbols <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_ids, columns = "SYMBOL", keytype = "ENTREZID")



if(nrow(gene_symbols)>0){
  for(i in 1:length(gene_symbols$ENTREZID)){
    gene_track@range@elementMetadata@listData$id[gene_track@range@elementMetadata@listData$gene == gene_symbols$ENTREZID[i]] <- gene_symbols$SYMBOL[i]
    gene_track@range@elementMetadata@listData$symbol[gene_track@range@elementMetadata@listData$gene == gene_symbols$ENTREZID[i]] <- gene_symbols$SYMBOL[i]
    
    
  }
}

plotTracks(gene_track,
           from =view_win[1],
           to=view_win[2]   )

plotTracks(gene_track,
from = min( plot_df$pos[which(plot_df$study=="AD_Bellenguez_2022")]),
to=max( plot_df$pos[which(plot_df$study=="AD_Bellenguez_2022")]) )
