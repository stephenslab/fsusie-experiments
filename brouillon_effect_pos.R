


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
                             ylim =c( min( c(effect_s)),max(c(effect_s)  )),
                             col = "royalblue", pch = 16, cex = 1.2,
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
  
 
}



tt= do.call( rbind , df_list)
 
idl= which( pos > view_win[1] & pos < view_win[2])-1
total_overlay= OverlayTrack( trackList =plot_list[idl],
                             background.title = "white")


 


effect0=       rep(0,length(out$effect_estimate ))
group_cred= c( 0)
group_colors <- c("black"  )

group_lwd= c(2)
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

fsusie_ha_plot <- OverlayTrack(trackList = list( haQTL_pos0,
                                                OverlayTrack(trackList =plot_list[idl])
                                                ),
                               background.title = "white"
)


list_track=  list(     
  otAD,
  otCR1,
  otCR2  ,
  
  t_ha,
  fsusie_ha_plot ,
  gene_track,
  gtrack 
  
)

plotTracks(list_track,
           from =view_win[1],
           to=view_win[2])

