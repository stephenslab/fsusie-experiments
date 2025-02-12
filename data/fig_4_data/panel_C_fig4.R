rm(list=ls())

library(AnnotationHub)
library(org.Hs.eg.db)
library(GenomicRanges)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(fsusieR)
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



##" cas 1 ---- 
plot_df <- data$d[[1]]
haQTL_df <- data$d[[2]]
mQTL_df <- data$d[[3]]
gene_info <- data$d[[4]]
sumstat <- data$d[[5]]
pip_df <- data$d[[6]]


# Parameters
view_win <- c(4759843, 5000000)
text_size <- 20
gene <- c("ENSG00000029725", "ENSG00000161929", "ENSG00000108556")
custom_labeller <- function(x) {
  x %>% 
    gsub("DeJager_", "", ., fixed = TRUE) %>%  
    gsub("([_:,|-])", "\n", .)             
}
# Function to extract SNP position from a given notation


plot_df <- data$d[[1]]
haQTL_df <- data$d[[2]]
mQTL_df <- data$d[[3]]
gene_info <- data$d[[4]]
sumstat <- data$d[[5]]
pip_df <- data$d[[6]]


# Parameters
view_win <- c(4759843, 5000000)








chr =  paste("chr", 17, sep = "")


data_track =plot_df[ which(  plot_df$study == "DLPFC_DeJager_eQTL" ),]
data_track = data_track[which(    data_track$pos > view_win[1] &  data_track$pos <view_win[2] ),]

data_track_CS1 =plot_df [ which(  plot_df$study == "DLPFC_DeJager_eQTL" & plot_df$CS1 ),]  #%>% filter(study == "DLPFC_DeJager_eQTL", CS1)#"maroon"
data_track_CS4 =plot_df [ which(  plot_df$study == "DLPFC_DeJager_eQTL" & plot_df$CS4 ),]# plot_df %>% filter(CS4, study == "DLPFC_DeJager_eQTL") #"steelblue"



t1= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track$pos , end = data_track$pos )),
                data = matrix(data_track$z , nrow=1), genome = "hg19", 
                type = "p", col = "black",  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black") ) # Change title color to black




t2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track_CS1$pos , end = data_track_CS1$pos )),
                data = matrix(data_track_CS1$z , nrow=1), genome = "hg19", 
                ylim =c( min(data_track$z), max(data_track$z)),
                type = "p", col = "maroon", cex=1.5,  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black") ) # Change title color to black


t3= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track_CS4$pos , end = data_track_CS4$pos )),
                data = matrix(data_track_CS4$z , nrow=1), genome = "hg19", 
                ylim =c( min(data_track$z), max(data_track$z)),
                type = "p", col = "steelblue",
                cex=1.5,# Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black") ) # Change title color to black


otAD <- OverlayTrack(trackList=list(    t1, t3 ,t2 ))


plotTracks( otAD)

# second panle

data_ha =pip_df[which( pip_df$study =="ROSMAP_DLPFC_haQTL"&pip_df$cs_coverage_0.95_min_corr==2  ),]
#  pip_df %>% filter(study %in% c("ROSMAP_DLPFC_haQTL"), cs_coverage_0.95_min_corr == 2)
data_ha= data_ha[which(data_ha$pos> view_win[1] & data_ha$pos<view_win[2]),]
t_ha= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_ha$pos , end = data_ha$pos )),
                  data = matrix(data_ha$pip , nrow=1), genome = "hg19", 
                  ylim =c( 0, 0.5),
                  type = "p", col = "steelblue",
                  cex=1.5,# Use color column from df_plot
                  track.margin = 0.05, # Reduce margin between track and title
                  cex.title = 0.6,     # Reduce title size
                  cex.axis = 0.6,      # Reduce axis text size
                  col.axis = "black",  # Change axis color to black
                  col.title = "black") ) # Change title color to black


plotTracks( t_ha)


data_me =  pip_df[which( pip_df$study %in% c("ROSMAP_DLPFC_mQTL", "") & pip_df$cs_coverage_0.95 == 7),]
  #pip_df %>% filter(study %in% c("ROSMAP_DLPFC_mQTL", ""), cs_coverage_0.95 == 7)
data_me= data_me[which(data_me$pos> view_win[1] & data_me$pos<view_win[2]),]
t_me= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_me$pos , end = data_me$pos )),
                  data = matrix(data_me$pip , nrow=1), genome = "hg19", 
                  ylim =c( 0, 0.5),
                  type = "p", col = "maroon",
                  cex=1.5,# Use color column from df_plot
                  track.margin = 0.05, # Reduce margin between track and title
                  cex.title = 0.6,     # Reduce title size
                  cex.axis = 0.6,      # Reduce axis text size
                  col.axis = "black",  # Change axis color to black
                  col.title = "black") ) # Change title color to black
otme_ha <- OverlayTrack(trackList=list(    t_me, t_ha  ))
plotTracks(otme_ha, from =view_win[1], to = view_win[2])

#### ha ---- ----- 

res <- readRDS(paste0(path, "/data/fig_4_data/fsusie_object/ROSMAP_haQTL.chr17_1059843_6175034.fsusie_mixture_normal_top_pc_weights.rds"))

fsusie_obj_ha=res$`chr17:1059843-6175034`$ROSMAP_DLPFC_haQTL$fsusie_result
rm(res)
positions = fsusie_obj_ha$outing_grid


effect=  fsusie_obj_ha$fitted_func[[2]]


haQTL_track = DataTrack(range = GRanges(seqnames = chr,
                                        ranges = IRanges(start = positions ,
                                                         end = positions   + 1)),
                        data = effect , genome = "hg38",
                        type = "l", 
                        track.margin = 0.05 ,
                        col.axis = "black",col.title = "black",
                        fontface = "plain",background.title = "white",
                        fontface.title = 1)



effect=  fsusie_obj_ha$cred_band[[2]][1, ]




haQTL_trackcb1  = DataTrack(range = GRanges(seqnames = chr,
                                            ranges = IRanges(start = positions ,
                                                             end = positions + 1)),
                            data = effect , genome = "hg38",
                            type = "l", 
                            track.margin = 0.05 ,lty=2,
                            col.axis = "black",col.title = "black",
                            fontface = "plain",background.title = "white",
                            fontface.title = 1)


effect=  fsusie_obj_ha$cred_band[[2]][2, ]




haQTL_trackcb2  = DataTrack(range = GRanges(seqnames = chr,
                                            ranges = IRanges(start = positions ,
                                                             end = positions + 1)),
                            data = effect , genome = "hg38",
                            type = "l", 
                            track.margin = 0.05 ,lty=2,
                            col.axis = "black",col.title = "black",
                            fontface = "plain",background.title = "white",
                            fontface.title = 1)

fsusie_ha_plot <- OverlayTrack(trackList=list( haQTL_track,haQTL_trackcb1, haQTL_trackcb2 ))
plotTracks(fsusie_ha_plot , from =view_win[1], to = view_win[2])


#### meqtl -----

res <- readRDS(paste0(path, "/data/fig_4_data/fsusie_object/ROSMAP_mQTL.chr17_1059843_6175034.fsusie_mixture_normal_top_pc_weights.rds"))
fsusie_obj_me = res$`chr17:1059843-6175034`$ROSMAP_DLPFC_mQTL$fsusie_result 



positions = fsusie_obj_me$outing_grid

out=list(effect= fsusie_obj_me$fitted_func[[7]],
         cred_band=  fsusie_obj_me$cred_band[[7]]  )
effect= out$effect


meQTL_track = DataTrack(range = GRanges(seqnames = chr,
                                        ranges = IRanges(start = positions ,
                                                         end = positions   + 1)),
                        data = effect , genome = "hg38",
                        type = "l", 
                        track.margin = 0.05 ,
                        col.axis = "black",col.title = "black",
                        fontface = "plain",background.title = "white",
                        fontface.title = 1)



effect=  out$cred_band[1,]




meQTL_trackcb1  = DataTrack(range = GRanges(seqnames = chr,
                                            ranges = IRanges(start = positions ,
                                                             end = positions + 1)),
                            data = effect , genome = "hg38",
                            type = "l", 
                            track.margin = 0.05 ,lty=2,
                            col.axis = "black",col.title = "black",
                            fontface = "plain",background.title = "white",
                            fontface.title = 1)


effect=    out$cred_band[2,]




meQTL_trackcb2  = DataTrack(range = GRanges(seqnames = chr,
                                            ranges = IRanges(start = positions ,
                                                             end = positions + 1)),
                            data = effect , genome = "hg38",
                            type = "l", 
                            track.margin = 0.05 ,lty=2,
                            col.axis = "black",col.title = "black",
                            fontface = "plain",background.title = "white",
                            fontface.title = 1)

fsusie_me_plot <- OverlayTrack(trackList=list( meQTL_track,meQTL_trackcb1, meQTL_trackcb2 ))
plotTracks(fsusie_me_plot  )



fsusie_me_plot <- OverlayTrack(trackList=list( haQTL_track,haQTL_trackcb1, haQTL_trackcb2,
                                               meQTL_track,meQTL_trackcb1, meQTL_trackcb2 ))
plotTracks(fsusie_me_plot , from =view_win[1], to = view_win[2])


fsusie_me_plot <- OverlayTrack(trackList=list( meQTL_track,meQTL_trackcb1, meQTL_trackcb2  ))
plotTracks(fsusie_me_plot , from =view_win[1], to = view_win[2])

 

genome_track <- GenomeAxisTrack(col.axis = "black",col.title = "black")

# Create a "gene region" track.
gene_track <- GeneRegionTrack(txdb,genome = "hg38",chromosome = chr,
                              pos0 = view_win[1],pos1 = view_win[2],name = "",
                              showId = TRUE,geneSymbol = TRUE,
                              col.axis = "black",col.title = "black",
                              transcriptAnnotation = "symbol",
                              rotation.title = 0 ,
                              col = "salmon",fill = "salmon",
                              background.title = "white")

plotTracks(gene_track)
tracks <- c(otAD,
  otme_ha,
  fsusie_me_plot,
  gene_track
)
plotTracks(gene_track , from =view_win[1], to = view_win[2])

