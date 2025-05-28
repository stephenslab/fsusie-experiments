source("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/scripts_plot/cases_study/get_gene_annotations.R", echo=TRUE)
source("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/scripts_plot/cases_study/interpolate_effect_estimates.R", echo=TRUE)

path="C:/Document/Serieux/Travail/Data_analysis_and_papers/GTEX_analysis_Fsusie"
source(paste0(path,"/code/plot_all_effect_log.R"))
source(paste0(path,"/code/plot_log.R"))
library(fsusieR)
library(ggplot2)

cs_colors <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99",
               "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#b15928",
               "#8dd3c7","#bebada","#fb8072","#80b1d3","#fdb462",
               "#b3de69","#fccde5","#bc80bd")

p0="C:/Document/Serieux/Travail/Data_analysis_and_papers/GTEX_analysis_Fsusie_res/result_log_local_scale"
load(paste0(p0,"/ENSG00000099194.csv.gz.RData"))
 
zoomin_region=c(out$locus[1],  out$locus[2])/1e6
ids <- names(out$res$pip)
pdat <- data.frame(id  = as.character(NA),
                   pos =as.numeric(sapply(strsplit(ids, "_"), function(x) x[2])),
                   pip = out$res$pip,
                   cs  = as.character(NA),
                   stringsAsFactors = FALSE)
pdat <- transform(pdat,pos = as.numeric(pos)/1e6)
rownames(pdat) <- ids
zoomout_region <- range(pdat$pos)
n <- length(out$res$cs)
for (i in 1:n) {
  snps <- names(out$res$cs[[i]])
  cs_label <- sprintf("CS %d, %d SNPs",i,length(snps))
  pdat[snps,"cs"] <- i # cs_label
  j <- snps[which.max(pdat[snps,"pip"])]
  pdat[j,"id"] <- sprintf("%s (CS %d, %d SNPs)",j,i,length(snps))
}


 



rownames(pdat) <- ids
zoomout_region <- range(pdat$pos)



n <- length(obj_plot$fsusie_obj_me$sets$cs)
for (i in 1:n) {
  snps <- names(obj_plot$fsusie_obj_me$sets$cs[[i]])
  cs_label <- sprintf("CS %d, %d SNPs",i,length(snps))
  pdat[snps,"cs"] <- i # cs_label
  j <- snps[which.max(pdat[snps,"pip"])]
  pdat[j,"id"] <- sprintf("%s (CS %d, %d SNPs)",j,i,length(snps))
}
pdat <- transform(pdat,cs = factor(cs))
i <- which(pdat$pip >= 0.01)
j <- which(pdat$pip < 0.01)
j <- sample(j,1000)
i <- sort(c(i,j))
pdat <- pdat[i,]
p1 <- ggplot(pdat, aes(x = pos, y = pip, label = id)) +
  geom_rect(data = data.frame(xmin = zoomin_region[1],
                              xmax = zoomin_region[2],
                              ymin = -0.05,
                              ymax = 1.1),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightblue", alpha = 0.4, inherit.aes = FALSE) +
  geom_point(color = "black", size = 0.5) +
  geom_point(data = subset(pdat, !is.na(cs)), aes(color = cs),
             shape = 1, size = 1.75) +
  geom_text_repel(color = "midnightblue", size = 2.25,
                  min.segment.length = 0, max.overlaps = Inf,
                  segment.color = "midnightblue") +
  scale_x_continuous(breaks = seq(100.1, 100.5, 0.1)) +
  scale_y_continuous(limits = c(-0.05, 1.2), breaks = c(0, 0.5, 1)) +
  scale_color_manual(values = cs_colors, na.value = "darkgray") +
  guides(color = "none") +
  labs(x = "base-pair position on chromosome 20 (Mb)", y = "PIP",
       title = "SCD") +
  theme_cowplot(font_size = 8) +
  theme(
    axis.line.x = element_blank(),
    axis.ticks = element_line(color = "black") 
    
  )
p1
