library(data.table)


Id_map = fread("C:/Document/Serieux/Travail/Data_analysis_and_papers/xQTL_data/H3K9ac_count/sampleSheetAfterQc.csv")


meta_data= fread("C:/Document/Serieux/Travail/Data_analysis_and_papers/xQTL_data/H3K9ac_count/H3K9acDomains.csv" )


tt= fread("C:/Document/Serieux/Travail/Data_analysis_and_papers/xQTL_data/H3K9ac_count/H3K9acCounts.txt" )
tt= t(tt)


view_win <- c(51200000, 51480000)
chr=12

idx = which( meta_data$chr==12 &  meta_data$start> view_win[1] & meta_data$start< view_win[2])

dim(tt)
peak_names= (tt[1,])
peak_counts= apply( tt[ -1,idx],2,as.numeric)
width=   meta_data[idx,]$ width 

l_name=rep(NA, nrow(peak_counts))
for (i in 2:  nrow(tt)){
  
   if ( length(Id_map$SampleID[which(Id_map$ProjID == row.names(tt)[i])] )>0){
      l_name[i-1]= Id_map$SampleID[which(Id_map$ProjID == row.names(tt)[i])] 
   }
   
}

rownames(peak_counts)=l_name
peak_counts= peak_counts [-which (is.na(rownames(peak_counts))),]

dim(peak_counts)


res <- readRDS("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/fsusie_object/ROSMAP_haQTL.chr12_50815042_54677408.fsusie_mixture_normal_top_pc_weights.rds")


res_ha <- readRDS("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/fsusie_object/raw_data/ROSMAP_haQTL.chr12_50815042_54677408.fsusie_mixture_normal_top_pc_weights.input_data.rds")



Y= as.data.frame(res_ha$residual_Y)


X=as.data.frame(res_ha$residual_X)

X=X[order(row.names(X)),]
peak_counts=peak_counts[order(rownames(peak_counts)),]

share_name= intersect(rownames(peak_counts),row.names(X))

X= X[which(row.names(X) %in% share_name),]

peak_counts=peak_counts[which(row.names(peak_counts) %in% share_name),]





lead_SNP= matrix(X[,2716], ncol=1)


col_names <-        colnames(as.data.frame(res_ha$X_data)) 

pos_SNP_HA <-  as.numeric(gsub("chr[0-9XY]+\\.([0-9]+)\\..*", "\\1", col_names))

pos = as.data.frame(res_ha$Y_coordinates) #use start
pos= pos$start


map_data <- fsusieR:::remap_data(Y=Y,
                                 pos=pos,
                                 
                                 max_scale=10)

outing_grid <- map_data$outing_grid
Y_w= map_data$Y






hist(lead_SNP, nclass=100)

binarized_SNP= ifelse(lead_SNP>0,1,0)

#51204562 -51204865

 



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
 





idx_SNP_0= which(binarized_SNP==0)

idx_SNP_1= which(binarized_SNP==1)

pos <-  0.5*meta_data$width[idx]  +  meta_data$start[idx]
height <- apply(peak_counts[idx_SNP_0,],2,mean ) # Random heights between 0.5 and 2
width <- meta_data$width[idx]   # Random widths between 0.5 and 2
id2= which( width>1000)
pos= pos[-id2]
height= height[-id2]
width= width[-id2]/2


x_range <- seq(min(pos) - 3 * max(width), max(pos) + 3 * max(width), length.out = 1000)

plot(x_range, rep(0, length(x_range)), type = "n", 
     xlab = "X", ylab = "Density", 
     main = " ", ylim = c(0, max(height)))

for (i in seq_along(pos)) {
  y <- height[i] * exp(-((x_range - pos[i])^2) / (2 * (width[i]^2)))
  lines(x_range, y, 
        #col = colors[i],
        lwd = 1)
}




pos <-  0.5*meta_data$width[idx]  +  meta_data$start[idx]
height <- apply(peak_counts[idx_SNP_1,],2,mean ) # Random heights between 0.5 and 2
width <- meta_data$width[idx]   # Random widths between 0.5 and 2
id2= which( width>1000)
pos= pos[-id2]
height= height[-id2]
width= width[-id2]/2


for (i in seq_along(pos)) {
  y <- height[i] * exp(-((x_range - pos[i])^2) / (2 * (width[i]^2)))
  lines(x_range, y, 
         col ="blue3",
        lwd = 1)
}










library(ggplot2)
library(dplyr)

# Compute density curves
compute_curves <- function(idx_SNP, color) {
  pos <- 0.5 * meta_data$width[idx] + meta_data$start[idx]
  height <- apply(peak_counts[idx_SNP, ], 2, mean)
  width <- meta_data$width[idx] / 2
  
  # Filter out very wide peaks
  valid_idx <- which(width <= 1000)
  pos <- pos[valid_idx]
  height <- height[valid_idx]
  width <- width[valid_idx]
  
  x_range <- seq(min(pos) - 3 * max(width), max(pos) + 3 * max(width), length.out = 1000)
  
  data_list <- list()
  
  for (i in seq_along(pos)) {
    y_values <- height[i] * exp(-((x_range - pos[i])^2) / (2 * (width[i]^2)))
    data_list[[i]] <- data.frame(x = x_range, y = y_values, group = color)
  }
  
  return(do.call(rbind, data_list))
}

# Generate data for both SNP groups
data_snp0 <- compute_curves(idx_SNP_0, "SNP_0")
data_snp1 <- compute_curves(idx_SNP_1, "SNP_1")

# Combine data
plot_data <- bind_rows(data_snp0, data_snp1)

# Plot using ggplot2 with increased transparency
ggplot(plot_data, aes(x = x, y = y, color = group, fill = group,, alpha = 0.15)) +
  geom_line(size = 1, alpha = 0.15) +
  geom_ribbon(aes(ymin = 0, ymax = y), alpha = 0.15) + # Increased transparency
  scale_color_manual(values = c("SNP_0" = "red3", "SNP_1" = "green")) +
  scale_fill_manual(values = c("SNP_0" = "red3", "SNP_1" = "green")) +
  labs(x = "X", y = "Density", title = "Density Curves for SNP Groups") +
  theme_minimal()






plot_data[,2]= log1p(plot_data[,2])


count_0= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = data_snp0$x , end =data_snp0$x )),
                     data = matrix(log1p(data_snp0$y ), nrow=1), genome = "hg38", 
                     ylim=c(0,  max(plot_data[,2]) ) ,
                     type = "h", col = "red3",
                     pch=25,
                     fill= "red3",
                     cex=1.5,# Use color column from df_plot
                     track.margin = 0.05, # Reduce margin between track and title
                     cex.title = 0.6,     # Reduce title size
                     cex.axis = 0.6,      # Reduce axis text size
                     col.axis = "black",  # Change axis color to black
                     col.title = "black",rotation.title = 90,cex.title = cex,
                     background.title = "white",name="observed") )
count_1= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start =data_snp1$x , end = data_snp1$x)),
                     data = matrix(log1p( data_snp1$y ), nrow=1), genome = "hg38", 
                     ylim=c(0,  max(plot_data[,2]) ),
                     type = "h", col = "green3" ,
                     pch=25,
                     fill= "green3" ,
                     cex=1.5,# Use color column from df_plot
                     track.margin = 0.05, # Reduce margin between track and title
                     cex.title = 0.6,     # Reduce title size
                     cex.axis = 0.6,      # Reduce axis text size
                     col.axis = "black",  # Change axis color to black
                     col.title = "black",rotation.title = 90,cex.title = cex,
                     background.title = "white",name="observed") )


count_track=OverlayTrack( list(count_0,
                               count_1),
                          ,
                          background.title = "white")
plotTracks(
  count_track)

list_track=  list( otAD,
                   otGALNT6,
                   otSLC4A8 ,
                   t_me,t_ha,
                   
                   fsusie_me_plot ,
                   fsusie_ha_plot,
                   count_track,
                   gene_track,
                   gtrack
)

plotTracks(list_track,
           from =view_win[1],
           to=51480000,  
           frame = TRUE,
           
           cex.main=1.2, cex.title = 1.)

#view_win <- c(5.12e7, 5.16e7) 



folder_path=  paste0(getwd(),
                     "/plot/GALNT6/"
)
file_path <- file.path(folder_path, "GALNT6.pdf")
pdf(file_path, width = 8.27, height = 11.69)  # A4 in inches



grid.text(
  "rs3782473",
  x = 0.7,
  y =0.45,
  gp = gpar(col = "black", fontsize = 10)
)
dev.off() 
