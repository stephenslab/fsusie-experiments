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
                rotation.title = 90,
                background.title = "white",name="AD") ) # Change title color to black




t2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track_CS1$pos , end = data_track_CS1$pos )),
                data = matrix(data_track_CS1$z , nrow=1), genome = "hg19", 
                ylim =c( min(data_track$z), max(data_track$z)),
                type = "p", col = "royalblue", cex=1.5,  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
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
                rotation.title = 90,
                background.title = "white",name="GALNT6") ) # Change title color to black




t2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track_CS1$pos , end = data_track_CS1$pos )),
                data = matrix(data_track_CS1$z , nrow=1), genome = "hg19", 
                ylim =c( min(data_track$z), max(data_track$z)),
                type = "p", col = "royalblue", cex=1.5,  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
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
                rotation.title = 90,
                background.title = "white",name="SLC4A8") ) # Change title color to black




t2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track_CS1$pos , end = data_track_CS1$pos )),
                data = matrix(data_track_CS1$z , nrow=1), genome = "hg19", 
                ylim =c( min(data_track$z), max(data_track$z)),
                type = "p", col = "royalblue", cex=1.5,  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
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
                  ylim =c( 0.5, 1.1),
                  type = "p", col = "royalblue",
                  cex=1.5,# Use color column from df_plot
                  track.margin = 0.05, # Reduce margin between track and title
                  cex.title = 0.6,     # Reduce title size
                  cex.axis = 0.6,      # Reduce axis text size
                  col.axis = "black",  # Change axis color to black
                  col.title = "black",rotation.title = 90,cex.title = cex,
                  background.title = "white",name="PIP \n H3k4a9ac") ) # Change title color to black



data_me =  QTL_data %>%
  filter(variant_id == "chr12:51362485:T:C", study == "MSBB_mQTL") %>%
  mutate(study = "dmr-QTL") %>%
  separate(col = variant_id, into = c("chrom", "pos"), remove = FALSE) %>%
  mutate(pos = as.numeric(pos))

t_me= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_me$pos , end = data_me$pos )),
                  data = matrix(data_me$pip , nrow=1), genome = "hg38", 
                  ylim =c( 0.5, 1.1) ,
                  type = "p", col = "maroon",
                  cex=1.5,# Use color column from df_plot
                  track.margin = 0.05, # Reduce margin between track and title
                  cex.title = 0.6,     # Reduce title size
                  cex.axis = 0.6,      # Reduce axis text size
                  col.axis = "black",  # Change axis color to black
                  col.title = "black",rotation.title = 90,cex.title = cex,
                  background.title = "white",name="PIP \n DNAm") ) # Change title color to black






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

res <- readRDS("D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/fsusie_object/ROSMAP_haQTL.chr12_50815042_54677408.fsusie_mixture_normal_top_pc_weights.rds")
fsusie_obj_ha=res$`chr12:50815042-54677408`$ROSMAP_DLPFC_haQTL$fsusie_result
 rm(res)
#positions = fsusie_obj_ha$outing_grid


#effect=  fsusie_obj_ha$fitted_func[[5]]


#haQTL_track = DataTrack(range = GRanges(seqnames = chr,
                                        #                                        ranges = IRanges(start = positions ,
#                                                      end = positions   + 1)),
##                        data = effect , genome = "hg38",
#                        type = "l",  col = "steelblue",
#                        track.margin = 0.05 ,
#                        col.axis = "black",col.title = "black",
#                        fontface = "plain",rotation.title = 90,cex.title = cex,
#                        background.title = "white",name="effect H3k9ac")



#effect=  fsusie_obj_ha$cred_band[[5]][1, ]




#haQTL_trackcb1  = DataTrack(range = GRanges(seqnames = chr,
#                                            ranges = IRanges(start = positions ,
#                                                             end = positions + 1)),
#                            data = effect , genome = "hg38",
#                            type = "l", col = "steelblue",
#                            track.margin = 0.05 ,lty=2,
#                            col.axis = "black",col.title = "black",
#                            fontface = "plain",rotation.title = 90,cex.title = cex,
#                            background.title = "white",name="effect H3k9ac")


#effect=  fsusie_obj_ha$cred_band[[5]][2, ]




#haQTL_trackcb2  = DataTrack(range = GRanges(seqnames = chr,
#                                            ranges = IRanges(start = positions ,
#                                                             end = positions + 1)),
#                            data = effect , genome = "hg38",
#                            type = "l", col = "steelblue",
#                            track.margin = 0.05 ,lty=2,
#                            col.axis = "black",col.title = "black",
#                            fontface = "plain",rotation.title = 90,cex.title = cex,
#                            background.title = "white",name="effect H3k9ac")


#fsusie_ha_plot <- OverlayTrack(trackList=list( haQTL_track,haQTL_trackcb1, haQTL_trackcb2 ), background.title = "white")
#plotTracks(fsusie_ha_plot , from =view_win[1], to = view_win[2]      )




 res_ha <- readRDS("D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/fsusie_object/raw_data/ROSMAP_haQTL.chr12_50815042_54677408.fsusie_mixture_normal_top_pc_weights.input_data.rds")
 Y= as.data.frame(res_ha$residual_Y)


 X=as.data.frame(res_ha$residual_X)
 pos = as.data.frame(res_ha$Y_coordinates) #use start
 pos= pos$start


map_data <- fsusieR:::remap_data(Y=Y,
                                 pos=pos,
                                 
                                 max_scale=10)

outing_grid <- map_data$outing_grid
Y_w= map_data$Y


out= fsusieR:::univariate_TI_regression(Y_w,X= matrix(X[,2716], ncol=1),alpha=0.01)

positions=outing_grid 

effect=rbind(out$effect_estimate,
             out$cred_band,
             rep(0,length(out$effect_estimate)))
group_cred= c(1:3,0)
group_colors <- c("black" ,"royalblue","steelblue","steelblue" )

group_lwd= c(1,2,1,1)
haQTL_track =   DataTrack(range = GRanges(seqnames = chr,
                                          ranges = IRanges(start = positions,
                                                           end = positions + 1)),
                          data = effect, genome = "hg38",
                          groups= group_cred,
                          lwd = group_lwd,
                          rotation.title = 90,
                          name ="Effect H3k9ac",
                          type = "l",
                          col = group_colors,
                          type = c(  "s" ),
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



haQTL_pos= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = pos  , end = pos )),
                  data = matrix(0*pos , nrow=1), genome = "hg38", 
                  ylim =c( min( c(effect)),max(c(effect)  )) ,
                  type = "p", 
                  col = "lightblue" ,
                  cex= 1,# Use color column from df_plot
                  track.margin = 0.05, # Reduce margin between track and title
                  cex.title = 0.6,     # Reduce title size
                  cex.axis = 0.6,      # Reduce axis text size
                  col.axis = "black",  # Change axis color to black
                  col.title = "black", 
                  background.title = "white",
                  name="PIP DNAm") ) # Change title color to black



plotTracks(haQTL_track, from =view_win[1], to = view_win[2]      )
fsusie_ha_plot <- OverlayTrack(trackList=list(   haQTL_pos, haQTL_track    ),
                background.title = "white", show.legend = c( FALSE, FALSE) )

 
  plotTracks(fsusie_ha_plot , from =view_win[1], to = view_win[2]    , show.legend = FALSE  )


#### meqtl -----

res <-  readRDS("D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/fsusie_object/MSBB_mQTL.chr12_47653211_53108261.fsusie_mixture_normal_top_pc_weights.rds")
  
  fsusie_obj_me = res$`chr12:47653211-53108261`$MSBB_mQTL$fsusie_result


rm(res)
#positions = fsusie_obj_me$outing_grid

#out=list(effect= fsusie_obj_me$fitted_func[[14]],
#         cred_band=  fsusie_obj_me$cred_band[[14]]  )

res_me <- readRDS("D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/fsusie_object/raw_data/MSBB_mQTL.chr12_47653211_53108261.fsusie_mixture_normal_top_pc_weights.input_data.rds")

Y= as.data.frame(res_me$residual_Y)


X=as.data.frame(res_me$residual_X)
pos = as.data.frame(res_me$Y_coordinates) #use start
pos= pos$start


map_data <- fsusieR:::remap_data(Y=Y,
                                 pos=pos,
                                 
                                 max_scale=10)

outing_grid <- map_data$outing_grid
Y_w= map_data$Y

fsusie_obj_me$cs[[14]]
out= fsusieR:::univariate_TI_regression(Y_w,X= matrix(X[,9812], ncol=1),alpha=0.01)
plot( out$effect_estimate)
lines(out$cred_band[1,])                                        
lines(out$cred_band[2,])  



positions=outing_grid 


positions=outing_grid 

effect=rbind(out$effect_estimate,
             out$cred_band,
             rep(0,length(out$effect_estimate)))
group_cred= c(1:3,0)
group_colors <- c("black" ,"royalblue","steelblue","steelblue" )
group_lwd= c(1,2,1,1)

meQTL_track =   DataTrack(range = GRanges(seqnames = chr,
                                          ranges = IRanges(start = positions,
                                                           end = positions + 1)),
                          data = effect, genome = "hg38",
                          groups= group_cred,
                          rotation.title = 90,
                          name ="Effect DNAm",
                          type = "l",
                          lwd= group_lwd,
                          col = group_colors,
                          type = c(  "s" ),
                          cex=1.5,# Use color column from df_plot
                          track.margin = 0.05, # Reduce margin between track and title
                          cex.title = 1.6,     # Reduce title size
                          cex.axis = 0.6,      # Reduce axis text size
                          col.axis = "black",  # Change axis color to black
                          col.title = "black", 
                          background.title = "white",
                          legend=FALSE
                          )
plotTracks(meQTL_track , type = c(  "s" ))
meQTL_pos= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = pos  , end = pos )),
                       data = matrix(0*pos , nrow=1), genome = "hg38", 
                       ylim =c( min( c(effect)),max(c(effect))) ,
                       type = "p", 
                        
                       pch=20,
                       col = "royalblue" ,
                       cex=  1,# Use color column from df_plot
                       track.margin = 0.05, # Reduce margin between track and title
                       cex.title = 0.6,     # Reduce title size
                       cex.axis = 0.6,      # Reduce axis text size
                       col.axis = "black",  # Change axis color to black
                       col.title = "black",
                       background.title = "white" ,name="PIP DNAm") ) # Change title color to black



fsusie_me_plot <- OverlayTrack(trackList=list(  meQTL_pos, meQTL_track  ),
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
           to=max( plot_df$pos[which(plot_df$study=="Oli_mega_eQTL")]) )

plotTracks(list_track,
           from =view_win[1],
           to=view_win[2])

 

 

plotTracks(list_track,
           from =view_win[1],
           to=52000000 , frame = TRUE)
                000000
 
library(biomaRt)
library(GenomicRanges)

# Define the genomic region
chr <- "chr12"
start_pos <-  min( plot_df$pos[which(plot_df$study=="AD_Bellenguez_2022")])
end_pos <- max( plot_df$pos[which(plot_df$study=="Oli_mega_eQTL")]) 

# Use biomaRt to fetch gene annotations from Ensembl
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# Get gene and transcript information
genes <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", 
                 "strand", "ensembl_gene_id", "ensembl_transcript_id", "external_gene_name"),
  filters = c("chromosome_name", "start", "end"),
  values = list("12", start_pos, end_pos),
  mart = mart
)

# Get exon-level information
exons <- getBM(
  attributes = c("chromosome_name", "exon_chrom_start", "exon_chrom_end",
                 "strand", "ensembl_gene_id", "ensembl_transcript_id", 
                 "ensembl_exon_id", "external_gene_name"),
  filters = c("chromosome_name", "start", "end"),
  values = list("12", start_pos, end_pos),
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

  


list_track=  list( otAD,
                   otGALNT6,
                   otSLC4A8 ,
                    t_me,t_ha,
                   
                   fsusie_me_plot ,
                   fsusie_ha_plot,
                   gene_track
)

#view_win <- c(5.12e7, 5.16e7) 
plotTracks(list_track,
           from =view_win[1],
           to=view_win[2])



plotTracks(list_track,
           from =view_win[1],
           to=view_win[2],  frame = TRUE)
plotTracks(list_track,
           from = min( plot_df$pos[which(plot_df$study=="AD_Bellenguez_2022")]),
           to=51875000 ,
           frame = TRUE,
         
           sizes = c(1, 1,1, 0.5, 0.5,1,1,1))

list_track=  list( otAD,
                   otGALNT6,
                   otSLC4A8 ,
                   #t_me,t_ha,
                   
                   fsusie_me_plot ,
                   fsusie_ha_plot,
                   gene_track
)

plotTracks(list_track,
           from = min( plot_df$pos[which(plot_df$study=="AD_Bellenguez_2022")]),
           to=51875000 ,
           frame = TRUE,
           
           sizes = c(0.5, 0.5,0.5,  1,1,0.75))

 