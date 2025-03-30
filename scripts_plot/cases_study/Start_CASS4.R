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

cex=1
path= getwd()

AD_GWAS <- readRDS("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/GWAS_sumstat.chr20_55439357_57610823.rds")[[2]]
qTLdata= fread(paste0("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/eQTL.chr20_55439357_57610823.tsv.tsv"))

chr=20
view_win=c(min(AD_GWAS$pos),max(AD_GWAS$pos))
pdat1= AD_GWAS
dim(pdat1)
#pdat1$X.log10.P. 



t1= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start =pdat1$pos ,
                                                                 end = pdat1$pos )),
                data = matrix( (pdat1$`-log10(P)`) , nrow=1), genome = "hg19",
                ylim =c( min( (pdat1$`-log10(P)`)), max( (pdat1$`-log10(P)`))+0.2),
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
                data = matrix( (pdat1CS$`-log10(P)`) , nrow=1), genome = "hg19", 
                ylim =c( min( (pdat1$`-log10(P)`)), max( (pdat1$`-log10(P)`))+0.2),
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
otAD <- OverlayTrack(trackList=list(    t1,  t2  ),
                     background.title = "white")

plotTracks( otAD , from= view_win[1],
            to= view_win[2])









## CASS4 panel ----- 




pdat2 = qTLdata
z_to_neglog10p <- function(z) {
  z <- abs(z)
  log_p <- pnorm(z, lower.tail = FALSE, log.p = TRUE)  # log(1 - Φ(z))
  log10_p <- log_p / log(10)
  -log10(2) - log10_p  # two-sided
} 
pdat2$X.log10.P. <- z_to_neglog10p(pdat2$z)


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
                background.title = "white",name="CASS4") ) # Change title color to black

plotTracks(t1, from= view_win[1],
           to= view_win[2])

pdat2CS = pdat2[which(! is.na( pdat2$cs_coverage_0.95)),]
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




#### reprendre ici  ------ 





# pip plot  -----










#### meqtl -----



res <- readRDS("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/ROSMAP_mQTL.chr20_53859688_57519449.fsusie_mixture_normal_top_pc_weights.rds")
fsusie_obj_me = res$`chr20:53859688-57519449`$ROSMAP_DLPFC_mQTL$fsusie_result
rm(res) 
res_me <- readRDS("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/ROSMAP_mQTL.chr20_53859688_57519449.fsusie_mixture_normal_top_pc_weights.input_data.rds")## to work from here
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
# we actually care about CS 13
#plot_susiF(fsusie_obj_me)
fsusie_obj_me$cs[[13]]

56408019 -8.05 chr20_55439357_57610823 TRUE  FALSE AD_Bell… 8.42e-16        15.1
2 chr20:56412160:G:C   20    56412160 -8.09 chr20_55439357_57610823 TRUE  FALSE AD_Bell… 6.04e-16        15.2
3 chr20:56413016:C:T   20    56413016 -8.04 chr20_55439357_57610823 TRUE  FALSE AD_Bell… 8.90e-16        15.1
4 chr20:56414777:T:A   20    56414777 -8.01 chr20_55439357_57610823 TRUE  FALSE AD_Bell… 1.18e-15        14.9
5 chr20:56423488:A:G   20    56423488


## reprendre de la WW---- 
plotTracks(list_track)
Y= as.data.frame(res_me$residual_Y)


X=as.data.frame(res_me$residual_X)
pos = as.data.frame(res_me$Y_coordinates) #use start
pos= pos$start-1000000#weird


map_data <- fsusieR:::remap_data(Y=Y,
                                 pos=pos,
                                 
                                 max_scale=10)

outing_grid <- map_data$outing_grid
Y_w= map_data$Y


out= fsusieR:::univariate_smash_regression(Y_w,X= matrix(X[,12212  ], ncol=1),alpha=0.01)








plot( out$effect_estimate)
lines(out$cred_band[1,])                                        

lines(out$cred_band[2,])       


positions=outing_grid 

effect_s=rbind(out$effect_estimate,
               out$cred_band,
               rep(0,length(out$effect_estimate)))



chrom=6
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
                                    col = "black", pch = 16, cex =  .7,
                                    background.title = "white" )
  }
  
  
}



tt= do.call( rbind , df_list)

idl=  1:(length(pos) ) #which( pos > view_win[1] & pos < view_win[2])-1
total_overlay= OverlayTrack( trackList =plot_list ,
                             background.title = "white")










effect0=       rep(0,length(out$effect_estimate ))
group_cred= c( 0)
group_colors <- c("black"  )

group_lwd= c(1)
meQTL_pos0 =   DataTrack(range = GRanges(seqnames = chr,
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

fsusie_me_plot <- OverlayTrack(trackList = list(meQTL_pos0,
                                                OverlayTrack(trackList =plot_list )
),
background.title = "white"
)

plotTracks(fsusie_me_plot  )





list_track=  list( otAD,
                   oteqTL,
                   pip_overlay ,
                   fsusie_me_plot 
)
plotTracks(list_track)





