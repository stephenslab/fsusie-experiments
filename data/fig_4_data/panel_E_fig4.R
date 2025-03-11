rm(list=ls())

path= getwd()
data = readRDS(paste0(path , 
                      "/data/fig_4_data/Fig4_data.rds"))

view_win <- c(207317782, 207895513)

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



otAD <- OverlayTrack(trackList=list(    t1,  t2  ),
                     background.title = "white")


AD_cs =data_track_CS1
plotTracks( otAD )




#### CR1 panel ---- 


gene_name="CR1"
study="DLPFC_DeJager_eQTL" 
data_track =plot_df[ which(  plot_df$study ==study & plot_df$region==gene_name ),]
#data_track = data_track[which(    data_track$pos > view_win[1] &  data_track$pos <view_win[2]& plot_df$region==gene_name  ),]

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
#data_track = data_track[which(    data_track$pos > view_win[1] &  data_track$pos <view_win[2]& plot_df$region==gene_name  ),]

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

res <- readRDS("D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/fsusie_object/ROSMAP_haQTL.chr1_205117782_208795513.fsusie_mixture_normal_top_pc_weights.rds")
fsusie_obj_ha = res$`chr1:205117782-208795513`$ROSMAP_DLPFC_haQTL$fsusie_result
rm(res)



fsusie_obj_ha$cs[[5]]


res_ha <- readRDS("D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/fsusie_object/raw_data/ROSMAP_haQTL.chr1_205117782_208795513.fsusie_mixture_normal_top_pc_weights.input_data (1).rds")
Y= as.data.frame(res_ha$residual_Y)


X=as.data.frame(res_ha$residual_X)
pos = as.data.frame(res_ha$Y_coordinates) #use start
pos= pos$start


map_data <- fsusieR:::remap_data(Y=Y,
                                 pos=pos,
                                 
                                 max_scale=10)

outing_grid <- map_data$outing_grid
Y_w= map_data$Y


out= fsusieR:::univariate_TI_regression(Y_w,X= matrix(X[,10765 ], ncol=1),alpha=0.01)
plot( out$effect_estimate)
lines(out$cred_band[1,])                                        

lines(out$cred_band[2,])       


positions=outing_grid 

effect_s=rbind(out$effect_estimate,
               out$cred_band,
               rep(0,length(out$effect_estimate)))



chrom=1
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
                                    col = "black", pch = 16, cex =  .41,
                                    background.title = "white" )
  }
  
  
}



tt= do.call( rbind , df_list)

idl= which( pos > view_win[1] & pos < view_win[2])-1
total_overlay= OverlayTrack( trackList =plot_list[idl],
                             background.title = "white")





effect0=       rep(0,length(out$effect_estimate ))
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


haQTL_pos00 =   DataTrack(range = GRanges(seqnames = chr,
                                          ranges = IRanges(start = positions,
                                                           end = positions  )),
                          data = 0* positions, genome = "hg38",
                          groups= group_cred,
                          ylim =c( min( c(effect_s)),max(c(effect_s)  )) ,
                          lwd = group_lwd,
                          rotation.title = 90,
                          name ="effect H3k9ac",
                          type = c(  "l" ),
                          col = group_colors,
                          cex=0.41,
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


fsusie_ha_plot <- OverlayTrack(trackList = list( haQTL_pos0,
                                                 OverlayTrack(trackList =plot_list[idl])
),
background.title = "white"
)

plotTracks(fsusie_ha_plot  )
 


##### PIP plots  -----

data_ha =pip_df[which( pip_df$study =="ROSMAP_DLPFC_haQTL"&pip_df$cs_coverage_0.95_min_corr==5  ),]
tdf= plot_df[ which(  plot_df$study ==study & plot_df$region==gene_name ),]
tdf$z=tdf$z*0 

col_names <-        colnames(as.data.frame(res_ha$X_data)) 

pos_SNP_HA <-  as.numeric(gsub("chr[0-9XY]+\\.([0-9]+)\\..*", "\\1", col_names))
### ici -----

t_dat= fsusie_obj_ha$pip
t_dat[fsusie_obj_ha$cs[[5]]]=t_dat[fsusie_obj_ha$cs[[5]]]+0.05

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
#  pip_df %>% filter(study %in% c("ROSMAP_DLPFC_haQTL"), cs_coverage_0.95_min_corr == 2)
#data_ha= data_ha[which(data_ha$pos> view_win[1] & data_ha$pos<view_win[2]),]
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



 
#207577223 207624893 207629207


 
t_ha=OverlayTrack(trackList=list( t_0, t_ha0, t_ha1 ),
                  background.title = "white")
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
                   gene_track,
  gtrack 
                 
)

#view_win <- c(5.12e7, 5.16e7) 
#plotTracks(list_track,
#           from =view_win[1],
#           to=view_win[2])

#plotTracks(list_track,
#           from = min( plot_df$pos[which( 
#             plot_df$study=="DLPFC_DeJager_eQTL")]),
#           to=max( plot_df$pos[which(plot_df$study=="AD_Bellenguez_2022")]) ,
#           frame = TRUE 
#)
 
#plotTracks(list_track,
#           from =206320523 ,
#           to=max( plot_df$pos[which(plot_df$study=="AD_Bellenguez_2022")]) ,
#           frame = TRUE ,
#           sizes = c(0.5,0.5, 0.5,0.3,  1, 0.75)
#)
 

 

folder_path=  paste0(getwd(),
                     "/plot/CR1_CR2/"
)
file_path <- file.path(folder_path, "CR1_CR2.pdf")
pdf(file_path, width = 8.27, height = 11.69)  # A4 in inches


plotTracks(list_track,
           from =view_win[1],
           to=view_win[2]  ,
           frame = TRUE ,
           sizes = c(0.75,0.75, 0.75,0.35,  0.75, 0.75,0.3),
           #fontsize  = 15
           cex.main=1.2, cex.title = 1.
)

grid.text(
  "rs679515",
  x = 0.45,
  y =0.45,
  gp = gpar(col = "black", fontsize = 10)
)
grid.text(
  "rs10863418",
  x = 0.55,
  y =0.4 ,
  gp = gpar(col = "black", fontsize = 10)
)
grid.text(
  "rs4844610",
  x = 0.64,
  y =0.44,
  gp = gpar(col = "black", fontsize = 10)
)

dev.off()



####  global PIP plot ---- 

col_names <-        colnames(as.data.frame(res_ha$X_data)) 

pos_SNP_HA <-  as.numeric(gsub("chr[0-9XY]+\\.([0-9]+)\\..*", "\\1", col_names))




plot_colors <-          c("black" , "steelblue4", 
                          "green4", "deeppink1",
                          "#6A3D9A","royalblue",  
                          "darkturquoise", "green1",
                          "yellow4")
SNP_in_cs = c(unlist( fsusie_obj_ha$cs)) 

idx= which( fsusie_obj_ha$pip  >0.05  )
fsusie_obj_ha$pip[idx[ -which(idx %in% SNP_in_cs  )] ]=0.0
pos_SNP =  pos_SNP_HA
obj     = fsusie_obj_ha
point_size = 1.25
L <- obj$L
y <- obj$pip
font_size = 10
col_y <- rep(0,length(y))
for (l in 1:L) {
  col_y[which(1:length(y) %in% obj$cs[[l]])] <- l
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

P_pip_ha



break

AD_cs
CR1_cs
CR2_cs
HA_cs

intersect(AD_cs$pos,
          HA_cs$pos )
#207577223 207624893 207629207
#rs679515  rs10863418  rs4844610

intersect(AD_cs$pos,
          CR1_cs$pos )
#207577223 207629207
intersect(AD_cs$pos,
                     CR2_cs$pos )
#207510847 207577223


 
intersect(CR1_cs$pos,
          HA_cs$pos )
#207573951 207577223 207611623 207612944 207613197 207613483 207629207
intersect(CR2_cs$pos,
             HA_cs$pos )
# 207573951 207577223 207598421 207611623 207612944 207613483 207621975 207625349 207625371 207626529

 
 
 intersect(CR2_cs$pos,
           CR1_cs$pos)
 #207564732 207573951 207577223 207611623 207612944 207613483