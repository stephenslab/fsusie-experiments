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


view_win <- c(5.12e7, 5.16e7)

### AD GWAS panel -----

chr =  paste("chr", 12, sep = "")



study="AD_Bellenguez_2022"
data_track = plot_df [ which(   plot_df $study == study   ),]
#data_track = data_track[which(    data_track$pos > view_win[1] &  data_track$pos <view_win[2]& plot_df$region==gene_name  ),]

data_track_CS1 =plot_df [ which(  plot_df$study == study & plot_df$CS1   ),]  #%>% filter(study == "DLPFC_DeJager_eQTL", CS1)#"maroon"



t1= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track$pos , end = data_track$pos )),
                data = matrix(data_track$`-log10(P)` , nrow=1), genome = "hg19", 
                type = "p", col = "black",  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
                background.title = "white",name="AD") ) # Change title color to black




t2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track_CS1$pos , end = data_track_CS1$pos )),
                data = matrix(data_track_CS1$`-log10(P)` , nrow=1), genome = "hg19", 
                ylim =c( min(data_track$`-log10(P)`), max(data_track$`-log10(P)`)),
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
                data = matrix(data_track$`-log10(P)` , nrow=1), genome = "hg19", 
                type = "p", col = "black",  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
                background.title = "white",name="GALNT6") ) # Change title color to black




t2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track_CS1$pos , end = data_track_CS1$pos )),
                data = matrix(data_track_CS1$`-log10(P)` , nrow=1), genome = "hg19", 
                ylim =c( min(data_track$`-log10(P)`), max(data_track$`-log10(P)`)),
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
                data = matrix(data_track$`-log10(P)` , nrow=1), genome = "hg19", 
                type = "p", col = "black",  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
                background.title = "white",name="SLC4A8") ) # Change title color to black




t2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track_CS1$pos , end = data_track_CS1$pos )),
                data = matrix(data_track_CS1$`-log10(P)` , nrow=1), genome = "hg19", 
                ylim =c( min(data_track$`-log10(P)`), max(data_track$`-log10(P)`)),
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
t_ha1= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_ha$pos , end = data_ha$pos )),
                  data = matrix(data_ha$pip , nrow=1), genome = "hg38", 
                  ylim =c( 0.5, 1.1),
                  type = "p", col = "royalblue",
                  cex=1.5,# Use color column from df_plot
                  track.margin = 0.05, # Reduce margin between track and title
                  cex.title = 0.6,     # Reduce title size
                  cex.axis = 0.6,      # Reduce axis text size
                  col.axis = "black",  # Change axis color to black
                  col.title = "black",rotation.title = 90,cex.title = cex,
                  background.title = "white",name="PIP \n H3k4a9ac") )
 
data_ha = pip_df %>%
  filter(study %in% c("ROSMAP_DLPFC_haQTL") )
#pip_df %>% filter(study %in% c("ROSMAP_DLPFC_mQTL", ""), cs_coverage_0.95 == 7) 
 



sub1= data_ha [which( data_ha$cs_coverage_0.95==4),]
t_ha2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start =sub1$pos , end =sub1$pos )),
                  data = matrix(sub1$pip , nrow=1), genome = "hg38", 
                  ylim =c( 0.0, 1.1),
                  type = "p", ,
                  col = c( "#6A3D9A"  ),
                  cex=1.5,# Use color column from df_plot
                  track.margin = 0.05, # Reduce margin between track and title
                  cex.title = 0.6,     # Reduce title size
                  cex.axis = 0.6,      # Reduce axis text size
                  col.axis = "black",  # Change axis color to black
                  col.title = "black",rotation.title = 90,cex.title = cex,
                  background.title = "white",name="PIP \n H3k4a9ac") ) # Change title color to black

plotTracks(t_ha1)

plotTracks(t_ha2)
 


t_ha= OverlayTrack(trackList = list(t_ha1,t_ha2),
                   background.title = "white")

plotTracks(t_ha )



data_me =  QTL_data %>%
  filter(variant_id == "chr12:51362485:T:C", study == "MSBB_mQTL") %>%
  mutate(study = "dmr-QTL") %>%
  separate(col = variant_id, into = c("chrom", "pos"), remove = FALSE) %>%
  mutate(pos = as.numeric(pos))

t_me2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_me$pos , end = data_me$pos )),
                  data = matrix(data_me$pip , nrow=1), genome = "hg38", 
                  ylim =c( 0.5, 1.1) ,
                  type = "p", col = "royalblue",
                  cex=1.5,# Use color column from df_plot
                  track.margin = 0.05, # Reduce margin between track and title
                  cex.title = 0.6,     # Reduce title size
                  cex.axis = 0.6,      # Reduce axis text size
                  col.axis = "black",  # Change axis color to black
                  col.title = "black",rotation.title = 90,cex.title = cex,
                  background.title = "white",name="PIP \n DNAm") ) # Change title color to black




data_me = pip_df %>%
  filter(study %in% c( "MSBB_mQTL") )  

sub= data_me [which(  data_me$pos > view_win[1] & data_me$pos < view_win[2]),]

sub1=data_me [which (data_me $cs_coverage_0.95==2),]
t_me1  = ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start =sub1$pos , end =sub1$pos )),
                    data = matrix(sub1$pip , nrow=1), genome = "hg38", 
                    ylim =c( 0,1 ) ,
                    type = "p", col = c(  "green4" ), 
                    cex=1.5,# Use color column from df_plot
                    track.margin = 0.05, # Reduce margin between track and title
                    cex.title = 0.6,     # Reduce title size
                    cex.axis = 0.6,      # Reduce axis text size
                    col.axis = "black",  # Change axis color to black
                    col.title = "black",rotation.title = 90,cex.title = cex,
                    background.title = "white",name="PIP \n DNAm") ) # Change title color to black
 


t_me  = OverlayTrack(trackList = list(t_me1,t_me2),
                     background.title = "white")

plotTracks(t_me )


list_track=  list( otAD,
                   otGALNT6,
                   otSLC4A8 ,
                   t_me,t_ha
                   
)
 

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
effect_s=rbind(out$effect_estimate,
               out$cred_band,
               rep(0,length(out$effect_estimate)))







chrom=12
plot_list=list()
df_list=list()
widthtick=1000
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
  ci_lower =   (1 - du/dupdi ) *effect_s[2, min(up)]+ (1 - di/dupdi ) *effect_s[2, max(low)]  
  ci_upper =  (1 - du/dupdi )* effect_s[3, min(up)]+ (1 - di/dupdi )* effect_s[3, max(low)]  
  df_list[[i-1]]= data.frame(effect=effect,
                             ci_lower=ci_lower,
                             ci_upper=ci_upper,
                             pos=pos[i])
  
  
  
  # DataTrack for effect size points (blue dots)
  dTrack_points <- DataTrack(start = pos[i ], end = pos[i ], genome = "hg38", chromosome = chrom,
                             data = effect, name = "Effect Size", type = "p",
                             ylim =c( -0.17,0.17  ) ,
                             col = "royalblue", pch = 16, cex =  .8,
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
                               col = "royalblue", pch = 16, cex = 1.2,
                               lwd=2,
                               ylim =c( -0.17,0.17  ) ,
                               
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
                       col = "royalblue", pch = 16, cex = 1.2,
                       ylim=c( -0.17,0.17  ) ,
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
                        col = "royalblue", pch = 16, cex = 1.2,
                        ylim =c( -0.17,0.17  ) 
  ) 
  
  
  tt = OverlayTrack(trackList = list( tick_low,tick_up,error_bar_track,  dTrack_points))
  #plotTracks(tt)
  
  plot_list[[i-1]] <- OverlayTrack(trackList = tt,
                                   background.title = "white")
  
  
}



tt= do.call( rbind , df_list)

idl= which( pos > view_win[1] & pos < view_win[2])-1
total_overlay1= OverlayTrack( trackList =plot_list[idl],
                              background.title = "white")



plotTracks(total_overlay1, from = view_win[1],
           to= view_win[2])

Y= as.data.frame(res_ha$residual_Y)


X=as.data.frame(res_ha$residual_X)
pos = as.data.frame(res_ha$Y_coordinates) #use start
pos= pos$start


map_data <- fsusieR:::remap_data(Y=Y,
                                 pos=pos,
                                 
                                 max_scale=10)

outing_grid <- map_data$outing_grid
Y_w= map_data$Y


out= fsusieR:::univariate_TI_regression(Y_w,X= matrix(X[,1991], ncol=1),alpha=0.01)

positions=outing_grid 
effect_s=rbind(out$effect_estimate,
               out$cred_band,
               rep(0,length(out$effect_estimate)))







chrom=12
plot_list=list()
df_list=list()
widthtick=1000
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
  ci_lower =   (1 - du/dupdi ) *effect_s[2, min(up)]+ (1 - di/dupdi ) *effect_s[2, max(low)]  
  ci_upper =  (1 - du/dupdi )* effect_s[3, min(up)]+ (1 - di/dupdi )* effect_s[3, max(low)]  
  df_list[[i-1]]= data.frame(effect=effect,
                             ci_lower=ci_lower,
                             ci_upper=ci_upper,
                             pos=pos[i])
  
  
  
  # DataTrack for effect size points (blue dots)
  dTrack_points <- DataTrack(start = pos[i ], end = pos[i ], genome = "hg38", chromosome = chrom,
                             data = effect, name = "Effect Size", type = "p",
                             ylim =  c( -0.17,0.17  ) ,
                             col = "#6A3D9A", pch = 16, cex =  .8,
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
                               col = "#6A3D9A", pch = 16, cex = 1.2,
                               lwd=2,
                               ylim = c( -0.17,0.17  ) ,
                               
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
                       
                       col = "#6A3D9A", pch = 16, cex = 1.2,
                       ylim = c( -0.17,0.17  ) ,
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
                        
                        col ="#6A3D9A", pch = 16, cex = 1.2,
                        ylim = c( -0.17,0.17  ) 
  ) 
  
  
  tt = OverlayTrack(trackList = list( tick_low,tick_up,error_bar_track,  dTrack_points))
  #plotTracks(tt)
  
  plot_list[[i-1]] <- OverlayTrack(trackList = tt,
                                   background.title = "white")
  
  
}



tt= do.call( rbind , df_list)

idl= which( pos > view_win[1] & pos < view_win[2])-1
total_overlay2= OverlayTrack( trackList =plot_list[idl],
                              background.title = "white")
plotTracks(total_overlay2, from = view_win[1],
           to= view_win[2])



effect0=       rep(0,length(out$effect_estimate ))
group_cred= c( 0)
group_colors <- c("black"  )

group_lwd= c(2)
haQTL_pos0 =   DataTrack(range = GRanges(seqnames = chr,
                                         ranges = IRanges(start = positions,
                                                          end = positions + 1)),
                         data = effect0, genome = "hg38",
                         groups= group_cred,
                         ylim =c( -0.17,0.17  ) ,
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

fsusie_ha_plot <- OverlayTrack(trackList = list( haQTL_pos0,
                                                 total_overlay2,
                                                 total_overlay1
),
background.title = "white"
)


plotTracks(fsusie_ha_plot )



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

effect_s=rbind(out$effect_estimate,
               out$cred_band,
               rep(0,length(out$effect_estimate)))







chrom=12
plot_listme=list()
df_list=list()
widthtick=1000
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
  ci_lower =   (1 - du/dupdi ) *effect_s[2, min(up)]+ (1 - di/dupdi ) *effect_s[2, max(low)]  
  ci_upper =  (1 - du/dupdi )* effect_s[3, min(up)]+ (1 - di/dupdi )* effect_s[3, max(low)]  
  df_list[[i-1]]= data.frame(effect=effect,
                             ci_lower=ci_lower,
                             ci_upper=ci_upper,
                             pos=pos[i])
  
  
  
  # DataTrack for effect size points (blue dots)
  dTrack_points <- DataTrack(start = pos[i ], end = pos[i ], genome = "hg38", chromosome = chrom,
                             data = effect, name = "Effect Size", type = "p",
                             ylim =c( min( c(effect_s)),max(c(effect_s)  )),
                             col = "royalblue", pch = 16, cex =  .8,
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
  
  plot_listme[[i-1]] <- OverlayTrack(trackList = tt,
                                     background.title = "white")
  
  
}



tt= do.call( rbind , df_list)

idlme= which( pos > view_win[1] & pos < view_win[2])-1
total_overlay1= OverlayTrack( trackList =plot_listme[idl],
                              background.title = "white")




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

effect_s=rbind(out$effect_estimate,
               out$cred_band,
               rep(0,length(out$effect_estimate)))







chrom=12
plot_listme=list()
df_list=list()
widthtick=1000
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
  ci_lower =   (1 - du/dupdi ) *effect_s[2, min(up)]+ (1 - di/dupdi ) *effect_s[2, max(low)]  
  ci_upper =  (1 - du/dupdi )* effect_s[3, min(up)]+ (1 - di/dupdi )* effect_s[3, max(low)]  
  df_list[[i-1]]= data.frame(effect=effect,
                             ci_lower=ci_lower,
                             ci_upper=ci_upper,
                             pos=pos[i])
  
  
  
  # DataTrack for effect size points (blue dots)
  dTrack_points <- DataTrack(start = pos[i ], end = pos[i ], genome = "hg38", chromosome = chrom,
                             data = effect, name = "Effect Size", type = "p",
                             ylim =c( min( c(effect_s)),max(c(effect_s)  )),
                             col = "royalblue", pch = 16, cex =  .8,
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
  
  plot_listme[[i-1]] <- OverlayTrack(trackList = tt,
                                     background.title = "white")
  
  
}



tt= do.call( rbind , df_list)

idlme= which( pos > view_win[1] & pos < view_win[2])-1
total_overlay2= OverlayTrack( trackList =plot_listme[idl],
                              background.title = "white")

