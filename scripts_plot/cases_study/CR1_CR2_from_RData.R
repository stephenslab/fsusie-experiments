 
rm(list=ls())

path= getwd()
load("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/plot/CR1_CR2/CR1_CR2_obj.RData")
# for writing the plots
folder_path=  paste0(getwd(),
                     "/plot/CR1_CR2/"
) 

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
bin_size=100
# Create a bin column
df$bin <- floor(df$obs_pos / bin_size)

# Load dplyr and aggregate
library(dplyr)

binned_df <- df %>%
  group_by(bin) %>%
  summarise(
    mean_func0 = mean(mean_func0, na.rm = TRUE),
    mean_func1 = mean(mean_func1, na.rm = TRUE),
    mean_func2 = mean(mean_func2, na.rm = TRUE),
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
