rm(list=ls())

library(AnnotationHub)
library(org.Hs.eg.db)
library(GenomicRanges)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

data = readRDS(  
                      "D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/Fig4_data.rds") 
extract_snp_position <- function(snp_string) {
  # Split the input string by ':'
  parts <- unlist(strsplit(snp_string, ":"))
  
  # Extract the position
  position <- as.numeric(parts[2])
  
  return(position)
}



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


data_track =plot_df %>% filter(study == "DLPFC_DeJager_eQTL" )
data_track = data_track[which(    data_track$pos > view_win[1] &  data_track$pos <view_win[2] ),]

data_track_CS1 =plot_df %>% filter(study == "DLPFC_DeJager_eQTL", CS1)#"maroon"
data_track_CS4 = plot_df %>% filter(CS4, study == "DLPFC_DeJager_eQTL") #"steelblue"


 
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


data_ha =  pip_df %>% filter(study %in% c("ROSMAP_DLPFC_haQTL"), cs_coverage_0.95_min_corr == 2)
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


data_me =  pip_df %>% filter(study %in% c("ROSMAP_DLPFC_mQTL", ""), cs_coverage_0.95 == 7)
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



data_effect_ha = haQTL_df %>% mutate(study = "haQTL effect") %>% filter(CS == 2)
library(dplyr)
library(tidyr)
df_reformatted <- data_effect_ha %>%
  pivot_longer(cols = c(start, end), names_to = "point", values_to = "x") %>%
  select(x, y = fun_plot)  # Rename fun_plot to y
library(dplyr)
library(tidyr)

# Reshape to ensure 'x' column doesn't conflict
df_reformatted <- data_effect_ha %>%
    pivot_longer(cols = c(start, end, x), names_to = "point", values_to = "x", names_repair = "minimal") 
plot( df_reformatted$fun_plot )

df_reformatted <-df_reformatted [  which(df_reformatted$point=="start"),]
df_reformatted= df_reformatted[order(df_reformatted$x),]



t_effect_ha= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = df_reformatted $x, end = df_reformatted $x )),
                  data = matrix(  df_reformatted $fun_plot , nrow=1), genome = "hg19", 
                   
                  type = "l", col = "maroon",
                  cex=1.5,# Use color column from df_plot
                  track.margin = 0.05, # Reduce margin between track and title
                  cex.title = 0.6,     # Reduce title size
                  cex.axis = 0.6,      # Reduce axis text size
                  col.axis = "black",  # Change axis color to black
                  col.title = "black") ) 


plotTracks(t_effect_ha, from =view_win[1], to = view_win[2])


 

data_effect_m = mQTL_df %>% mutate(study = "haQTL effect") %>% filter(CS == 7)


tt = list()
tt[[1]]=rbind( c( data_effect_ha$start[1], data_effect_ha$fun_plot[1]) ,
              c( data_effect_ha$end[1], data_effect_ha$fun_plot[1]))
for ( i in 2: (nrow(data_effect_m)-1 )){
  
  
  tt[[i]] =  rbind ( c( data_effect_ha$start[i], data_effect_ha$fun_plot[i]) ,
                    c( data_effect_ha$end[i], data_effect_ha$fun_plot[i+1]))
 
}

do.call( rbind,tt)

t_effect_m= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_effect_m$start, end =data_effect_m$end )),
                         data = matrix(data_effect_m$fun_plot , nrow=1), genome = "hg19", 
                         
                         type = "l", col = "maroon",
                         cex=1.5,# Use color column from df_plot
                         track.margin = 0.05, # Reduce margin between track and title
                         cex.title = 0.6,     # Reduce title size
                         cex.axis = 0.6,      # Reduce axis text size
                         col.axis = "black",  # Cmnge axis color to black
                         col.title = "black") ) 

plotTracks(t_effect_m , from =view_win[1], to = view_win[2])





data_me =       data = pip_df %>% filter(study %in% c("ROSMAP_DLPFC_mQTL", ""))


plot_df <- data$f[[1]]
haQTL_df <- data$f[[2]]
MSBB_df <- data$f[[3]]
gene_info <- data$f[[4]]
sumstat <- data$f[[5]]
pip_df <- data$f[[6]]
QTL_data <- data$f[[7]] 







### AD bellenguez ----- 
view_win <- c(5.12e7, 5.16e7)
data_track = plot_df[ which( plot_df$pos> view_win[1] & plot_df$pos < view_win[2]),  ] 

data_track = data_track[which( data_track$study=="AD_Bellenguez_2022"),]
chr =  paste("chr", 12, sep = "")
positions <-data_track$pos
t1= DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = positions, end = positions)),
          data = matrix(data_track$z , nrow=1), genome = "hg19", 
          type = "p", col = "black",  # Use color column from df_plot
          track.margin = 0.05, # Reduce margin between track and title
          cex.title = 0.6,     # Reduce title size
          cex.axis = 0.6,      # Reduce axis text size
          col.axis = "black",  # Change axis color to black
          col.title = "black")




view_win <- c(5.12e7, 5.16e7)
data_track = plot_df %>% filter(CS1)

data_track = data_track[which( data_track$study=="AD_Bellenguez_2022"),]
chr =  paste("chr", 12, sep = "")
positions <-data_track$pos
t2= DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = positions, end = positions)),
              data = matrix(data_track$z , nrow=1), genome = "hg19", 
              cex=2.5,
              type = "p", col = "steelblue",  # Use color column from df_plot
              track.margin = 0.05, # Reduce margin between track and title
              cex.title = 0.6,     # Reduce title size
              cex.axis = 0.6,      # Reduce axis text size
              col.axis = "black",  # Change axis color to black
              col.title = "black")


otAD <- OverlayTrack(trackList=list(t1, t2))


plotTracks( otAD)






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
extract_snp_position <- function(snp_string) {
  # Split the input string by ':'
  parts <- unlist(strsplit(snp_string, ":"))
  
  # Extract the position
  position <- as.numeric(parts[2])
  
  return(position)
}

# Example usage
snp_string <- "chr17:4899022:T:C"
position <- extract_snp_position(snp_string)
print(position)


 

















data_track =plot_df %>% filter(study == "DLPFC_DeJager_eQTL" )
data_track = data_track[which(    data_track$pos > view_win[1] &  data_track$pos <view_win[2] ),]

data_track_CS1 =plot_df %>% filter(study == "DLPFC_DeJager_eQTL", CS1)#"maroon"
data_track_CS4 = plot_df %>% filter(CS4, study == "DLPFC_DeJager_eQTL") #"steelblue"

groups= c(rep( "", nrow(data_track)),
          rep( "CS1", nrow(data_track_CS1)),
          rep( "CS4", nrow(data_track_CS4))
)
col_groups= c(rep ( "black", nrow(data_track) ),
              rep( "steelblue", nrow(data_track_CS1)),
              rep( "maroon", nrow(data_track_CS4))
              )

df_plot = rbind(data_track,
                data_track_CS1,data_track_CS4)

positions <-  do.call( c,
                       lapply( 1:nrow( df_plot),
                                  function( i) 
                                    extract_snp_position(df_plot$variant_alternate_id[i])
                               )
                       )

df_plot = rbind(data_track, data_track_CS1, data_track_CS4)

df_plot$color <- c(rep("black", nrow(data_track)),
                   rep("steelblue", nrow(data_track_CS1)),
                   rep("maroon", nrow(data_track_CS4)))
plotTracks( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = positions, end = positions)),
                      data = matrix(df_plot$z , nrow=1), genome = "hg19", 
                      type = "p", col = df_plot$color,  # Use color column from df_plot
                      track.margin = 0.05, # Reduce margin between track and title
                      cex.title = 0.6,     # Reduce title size
                      cex.axis = 0.6,      # Reduce axis text size
                      col.axis = "black",  # Change axis color to black
                      col.title = "black") ) # Change title color to black




plotTracks( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = positions, end =  positions)),
          data = matrix(df_plot$z , nrow=1) , genome = "hg19", 
          
          type = "p",col = col_groups  ,
          track.margin = 0.05, # Reduce margin between track and title
          cex.title = 0.6,     # Reduce title size
          cex.axis = 0.6,      # Reduce axis text size
          col.axis = "black",  # Change axis color to black
          col.title = "black") )# Change title color to black)











data_track = haQTL_df %>% mutate(study = "haQTL effect") %>% filter(CS == 2)
