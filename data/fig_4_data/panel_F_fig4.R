rm(list=ls())
library(ggplot2)
library(tidyr)
library(AnnotationHub)
library(org.Hs.eg.db)
library(GenomicRanges)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(fsusieR)
library(dplyr)
path= getwd()
data = readRDS(paste0(  path,"/data/fig_4_data/Fig4_data.rds") )
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


view_win <- c(4759843, 5300000) 

plot_df <- data$f[[1]]
haQTL_df <- data$f[[2]]
MSBB_df <- data$f[[3]]
gene_info <- data$f[[4]]
sumstat <- data$f[[5]]
pip_df <- data$f[[6]]
QTL_data <- data$f[[7]] 

 

### AD GWAS panel -----

chr =  paste("chr", 12, sep = "")

 

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
                col.title = "black",cex.title = cex,
                rotation.title = 0,
                background.title = "white",name="AD") ) # Change title color to black




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
                background.title = "white",name="AD") ) # Change title color to black



otAD <- OverlayTrack(trackList=list(    t1,  t2 ),
                     background.title = "white")



plotTracks( otAD )






#### GALNT6 panel ---- 
 

gene_name="GALNT6"
study="Oli_mega_eQTL" 
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
                background.title = "white",name="GALNT6") ) # Change title color to black




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
                background.title = "white",name="GALNT6") ) # Change title color to black
 

otGALNT6 <- OverlayTrack(trackList=list(    t1, t2 ),
                         background.title = "white")


plotTracks( otGALNT6 )










#### SLC4A8 panel ---- 

chr =  paste("chr", 12, sep = "")

gene_name="SLC4A8"
study="Oli_mega_eQTL" 
data_track =plot_df[ which(  plot_df$study ==study & plot_df$region==gene_name ),]
#data_track = data_track[which(    data_track$pos > view_win[1] &  data_track$pos <view_win[2]& plot_df$region==gene_name  ),]

data_track_CS1 =plot_df [ which(  plot_df$study == study & plot_df$CS1& plot_df$region==gene_name  ),]   


t1= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track$pos , end = data_track$pos )),
                data = matrix(data_track$z , nrow=1), genome = "hg19", 
                type = "p", col = "black",  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 0,
                background.title = "white",name="SLC4A8") ) # Change title color to black




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
                background.title = "white",name="SLC4A8"))  # Change title color to black



otSLC4A8 <- OverlayTrack(trackList=list(    t1,  t2 ),
                         background.title = "white")


plotTracks( otSLC4A8 )



####  PIP plot ------


data_ha = pip_df %>%
  filter(study %in% c("ROSMAP_DLPFC_haQTL"), cs_coverage_0.95 == 5)
#pip_df %>% filter(study %in% c("ROSMAP_DLPFC_mQTL", ""), cs_coverage_0.95 == 7) 
t_ha= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_ha$pos , end = data_ha$pos )),
                  data = matrix(data_ha$pip , nrow=1), genome = "hg38", 
                  ylim =c( 0, 1),
                  type = "p", col = "steelblue",
                  cex=1.5,# Use color column from df_plot
                  track.margin = 0.05, # Reduce margin between track and title
                  cex.title = 0.6,     # Reduce title size
                  cex.axis = 0.6,      # Reduce axis text size
                  col.axis = "black",  # Change axis color to black
                  col.title = "black",rotation.title = 0,cex.title = cex,
                  background.title = "white",name="PIP H3k4a9ac") ) # Change title color to black



data_me =  QTL_data %>%
  filter(variant_id == "chr12:51362485:T:C", study == "MSBB_mQTL") %>%
  mutate(study = "dmr-QTL") %>%
  separate(col = variant_id, into = c("chrom", "pos"), remove = FALSE) %>%
  mutate(pos = as.numeric(pos))

t_me= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_me$pos , end = data_me$pos )),
                  data = matrix(data_me$pip , nrow=1), genome = "hg38", 
                  ylim =c( 0, 1) ,
                  type = "p", col = "maroon",
                  cex=1.5,# Use color column from df_plot
                  track.margin = 0.05, # Reduce margin between track and title
                  cex.title = 0.6,     # Reduce title size
                  cex.axis = 0.6,      # Reduce axis text size
                  col.axis = "black",  # Change axis color to black
                  col.title = "black",rotation.title = 0,cex.title = cex,
                  background.title = "white",name="PIP DNAm") ) # Change title color to black






list_track=  list( otAD,
                   otGALNT6,
                   otSLC4A8 ,
                   t_me,t_ha
                   
)

view_win <- c(5.12e7, 5.16e7)
plotTracks(list_track,
           from = min( plot_df$pos[which(plot_df$study=="AD_Bellenguez_2022")]),
           to=max( plot_df$pos[which(plot_df$study=="AD_Bellenguez_2022")]) )

plotTracks(list_track,
           from =view_win[1],
           to=view_win[2])


#### ha ---- ----- 

res  <- readRDS("D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/fsusie_object/ROSMAP_haQTL.chr12_47653211_53108261.fsusie_mixture_normal_top_pc_weights.rds")

fsusie_obj_ha=res$`chr12:47653211-53108261`$ROSMAP_DLPFC_haQTL$fsusie_result
rm(res)
positions = fsusie_obj_ha$outing_grid


effect=  fsusie_obj_ha$fitted_func[[5]]


haQTL_track = DataTrack(range = GRanges(seqnames = chr,
                                        ranges = IRanges(start = positions ,
                                                         end = positions   + 1)),
                        data = effect , genome = "hg38",
                        type = "l",  col = "steelblue",
                        track.margin = 0.05 ,
                        col.axis = "black",col.title = "black",
                        fontface = "plain",rotation.title = 0,cex.title = cex,
                        background.title = "white",name="effect H3k9ac")



effect=  fsusie_obj_ha$cred_band[[5]][1, ]




haQTL_trackcb1  = DataTrack(range = GRanges(seqnames = chr,
                                            ranges = IRanges(start = positions ,
                                                             end = positions + 1)),
                            data = effect , genome = "hg38",
                            type = "l", col = "steelblue",
                            track.margin = 0.05 ,lty=2,
                            col.axis = "black",col.title = "black",
                            fontface = "plain",rotation.title = 0,cex.title = cex,
                            background.title = "white",name="effect H3k9ac")


effect=  fsusie_obj_ha$cred_band[[5]][2, ]




haQTL_trackcb2  = DataTrack(range = GRanges(seqnames = chr,
                                            ranges = IRanges(start = positions ,
                                                             end = positions + 1)),
                            data = effect , genome = "hg38",
                            type = "l", col = "steelblue",
                            track.margin = 0.05 ,lty=2,
                            col.axis = "black",col.title = "black",
                            fontface = "plain",rotation.title = 0,cex.title = cex,
                            background.title = "white",name="effect H3k9ac")


fsusie_ha_plot <- OverlayTrack(trackList=list( haQTL_track,haQTL_trackcb1, haQTL_trackcb2 ), background.title = "white")
plotTracks(fsusie_ha_plot , from =view_win[1], to = view_win[2]      )


#### meqtl -----

res <- readRDS("D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/fsusie_object/ROSMAP_mQTL.chr12_47653211_53108261.fsusie_mixture_normal_top_pc_weights (1).rds")
fsusie_obj_me = res$`chr12:47653211-53108261`$ROSMAP_DLPFC_mQTL$fsusie_result


rm(res)
positions = fsusie_obj_me$outing_grid

out=list(effect= fsusie_obj_me$fitted_func[[14]],
         cred_band=  fsusie_obj_me$cred_band[[14]]  )
effect= out$effect


meQTL_track = DataTrack(range = GRanges(seqnames = chr,
                                        ranges = IRanges(start = positions ,
                                                         end = positions   + 1)),
                        data = effect , genome = "hg38",
                        type = "l",  col = "maroon",
                        track.margin = 0.05 ,
                        col.axis = "black",col.title = "black",
                        fontface = "plain",rotation.title = 0,cex.title = cex,
                        background.title = "white",name="effect DNAm")



effect=  out$cred_band[1,]




meQTL_trackcb1  = DataTrack(range = GRanges(seqnames = chr,
                                            ranges = IRanges(start = positions ,
                                                             end = positions + 1)),
                            data = effect , genome = "hg38",
                            type = "l",  col = "maroon",
                            track.margin = 0.05 ,lty=2,
                            col.axis = "black",col.title = "black",
                            fontface = "plain",rotation.title = 0,cex.title = cex,
                            background.title = "white",name="effect DNAm")



effect=    out$cred_band[2,]




meQTL_trackcb2  = DataTrack(range = GRanges(seqnames = chr,
                                            ranges = IRanges(start = positions ,
                                                             end = positions + 1)),
                            data = effect , genome = "hg38",
                            type = "l",  col = "maroon",
                            track.margin = 0.05 ,lty=2,
                            col.axis = "black",col.title = "black",
                            fontface = "plain",rotation.title = 0,cex.title = cex,
                            background.title = "white",name="effect DNAm")


fsusie_me_plot <- OverlayTrack(trackList=list( meQTL_track,meQTL_trackcb1, meQTL_trackcb2 ),
                               background.title = "white")
plotTracks(fsusie_me_plot  )



#fsusie_me_plot <- OverlayTrack(trackList=list( haQTL_track,haQTL_trackcb1, haQTL_trackcb2,
#                                               meQTL_track,meQTL_trackcb1, meQTL_trackcb2 ),
#                               background.title = "white")
#plotTracks(fsusie_me_plot , from =view_win[1], to = view_win[2])
 




list_track=  list( otAD,
                   otGALNT6,
                   otSLC4A8 ,
                   t_me,t_ha,
                   
                   fsusie_me_plot ,
                   fsusie_ha_plot
)

view_win <- c(5.12e7, 5.16e7)
plotTracks(list_track,
           from = min( plot_df$pos[which(plot_df$study=="AD_Bellenguez_2022")]),
           to=max( plot_df$pos[which(plot_df$study=="AD_Bellenguez_2022")]) )

plotTracks(list_track,
           from =view_win[1],
           to=view_win[2])

plotTracks(list_track )

plotTracks(fsusie_me_plot,
           from =view_win[1],
           to=view_win[2])



view_win <- c(4759843, 5000000) 
plotTracks(list_track,
           from =view_win[1],
           to=view_win[2])


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

chr=paste0("chr",12)
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
                              background.title = "white")
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
           to=view_win[2])