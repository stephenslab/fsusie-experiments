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
path= getwd()


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


view_win <- c(56300000, 56700000) 
 
plot_df <- data$f[[1]]
haQTL_df <- data$f[[2]]
MSBB_df <- data$f[[3]]
gene_info <- data$f[[4]]
sumstat <- data$f[[5]]
pip_df <- data$f[[6]]
QTL_data <- data$f[[7]] 


chr =  paste("chr",20, sep = "")
 



### AD GWAS panel -----


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

plotTracks(t1)


t2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track_CS1$pos[1]  , end = data_track_CS1$pos[1]  )),
                data = matrix(data_track_CS1$`-log10(P)`[1] , nrow=1), genome = "hg19", 
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


t3= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_track_CS1$pos[-1] , end = data_track_CS1$pos[-1] )),
                data = matrix(data_track_CS1$`-log10(P)`[-1] , nrow=1), genome = "hg19", 
                ylim =c( min(data_track$`-log10(P)`), max(data_track$`-log10(P)`)+0.2),
                type = "p", col = "black", cex=1.5,
                fill=  "royalblue",
                pch=c(24),# Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
                background.title = "white",name="AD") ) 

tidx=which(data_track_CS1$pos==51354542)


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

otAD <- OverlayTrack(trackList=list(    t1,  t2,t3 ),
                     background.title = "white")



plotTracks( otAD )




res <- readRDS("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/CD2AP/ROSMAP_haQTL.chr20_53859688_57519449.fsusie_mixture_normal_top_pc_weights.input_data.rds")






##### effect plot  -----



res_ha <- readRDS("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/CD2AP/ROSMAP_haQTL.chr20_53859688_57519449.fsusie_mixture_normal_top_pc_weights.input_data.rds")


Y= as.data.frame(res_ha$residual_Y)


X=as.data.frame(res_ha$residual_X)
pos = as.data.frame(res_ha$Y_coordinates) #use start
pos= pos$start


map_data <- fsusieR:::remap_data(Y=Y,
                                 pos=pos,
                                 
                                 max_scale=10)

outing_grid <- map_data$outing_grid
Y_w= map_data$Y

#### change SNP index ------
out= fsusieR:::univariate_TI_regression(Y_w,X= matrix(X[,10765 ], ncol=1),alpha=0.01)
plot( out$effect_estimate)
lines(out$cred_band[1,])                                        

lines(out$cred_band[2,])       


positions=outing_grid 

effect_s=rbind(out$effect_estimate,
               out$cred_band,
               rep(0,length(out$effect_estimate)))



chrom=20
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



