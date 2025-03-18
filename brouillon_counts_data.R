library(data.table)


Id_map = fread("D:/Document/Serieux/Travail/Data_analysis_and_papers/xQTL_data/H3K9ac_count/sampleSheetAfterQc.csv")



tt= fread("D:/Document/Serieux/Travail/Data_analysis_and_papers/xQTL_data/H3K9ac_count/H3K9acCounts.txt" )
tt= t(tt)

dim(tt)
meta_data= fread("D:/Document/Serieux/Travail/Data_analysis_and_papers/xQTL_data/H3K9ac_count/H3K9acDomains.csv" )

res <- readRDS("D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/fsusie_object/ROSMAP_haQTL.chr12_50815042_54677408.fsusie_mixture_normal_top_pc_weights.rds")


res_ha <- readRDS("D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/fsusie_object/raw_data/ROSMAP_haQTL.chr12_50815042_54677408.fsusie_mixture_normal_top_pc_weights.input_data.rds")



Y= as.data.frame(res_ha$residual_Y)


X=as.data.frame(res_ha$residual_X)


col_names <-        colnames(as.data.frame(res_ha$X_data)) 

pos_SNP_HA <-  as.numeric(gsub("chr[0-9XY]+\\.([0-9]+)\\..*", "\\1", col_names))

pos = as.data.frame(res_ha$Y_coordinates) #use start
pos= pos$start


map_data <- fsusieR:::remap_data(Y=Y,
                                 pos=pos,
                                 
                                 max_scale=10)

outing_grid <- map_data$outing_grid
Y_w= map_data$Y


out= fsusieR:::univariate_TI_regression(Y_w,X= matrix(X[,2716], ncol=1),alpha=0.01)

meta_data


view_win <- c(5.12e7, 5.16e7)
chr=12

idx = which( meta_data$chr==12 &  meta_data$start> view_win[1] & meta_data$start< view_win[2])
lead_SNP= matrix(X[,2716], ncol=1)

hist(lead_SNP, nclass=100)

binarized_SNP= ifelse(lead_SNP>0,1,0)

51204562 -51204865
peak_names= (tt[1,])
 peak_counts= apply( tt[-1,idx],2,as.numeric)
width=   meta_data[idx,]$ width
hist( peak_counts[,34],nclass = 100)


plot(  , )
 


rownames(tt)



plot_multiple_bell_curves <- function(pos, height, width) {
  if (length(pos) != length(height) || length(pos) != length(width)) {
    stop("pos, height, and width must have the same length.")
  }
  
  x_range <- seq(min(pos) - 3 * max(width), max(pos) + 3 * max(width), length.out = 1000)
  
  # Create an empty plot
  plot(x_range, rep(0, length(x_range)), type = "n", 
       xlab = "X", ylab = "Density", 
       main = " ", ylim = c(0, max(height)))
  
  # Generate bell curves at different positions with different heights and widths
 # colors <- rainbow(length(pos))  # Different colors for each curve
  for (i in seq_along(pos)) {
    y <- height[i] * exp(-((x_range - pos[i])^2) / (2 * (width[i]^2)))
    lines(x_range, y, 
          #col = colors[i],
          lwd = 1)
  }
   
}
 


pos <-  0.5*meta_data$width[idx]  +  meta_data$start[idx]
height <- apply(peak_counts,2,mean ) # Random heights between 0.5 and 2
width <- meta_data$width[idx]   # Random widths between 0.5 and 2
id2= which( width>1000)
plot_multiple_bell_curves(pos[-id2], height[-id2], width[-id2]/2)

