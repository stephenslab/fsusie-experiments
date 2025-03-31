# Create a PIP "zoomout" plot for the CR1/CR2 locus.
library(data.table)
library(ggplot2)
library(ggrepel)
library(cowplot)
# Colors are from colorbrewer2.org
cs_colors <- c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e")
load("../outputs/CR1_CR2_obj.RData")
ids <- names(obj_plot$pip_fsusie_obj)
pdat <- data.frame(id  = as.character(NA),
                   pos = sapply(strsplit(ids,":",fixed = TRUE),"[[",2),
                   pip = obj_plot$pip_fsusie_obj,
                   cs  = as.character(NA),
                   stringsAsFactors = FALSE)
pdat <- transform(pdat,pos = as.numeric(pos)/1e6)
rownames(pdat) <- ids
n <- length(obj_plot$cs_fsusie_obj)
for (i in 1:n) {
  snps <- names(obj_plot$cs_fsusie_obj[[i]])
  cs_label <- sprintf("CS %d (%d SNPs)",i,length(snps))
  pdat[snps,"cs"] <- cs_label
  j <- snps[which.max(pdat[snps,"pip"])]
  pdat[j,"id"] <- j
}
pdat <- transform(pdat,cs = factor(cs))
p <- ggplot(pdat,aes(x = pos,y = pip,label = id)) +
  geom_point(color = "black",size = 0.5) +
  geom_point(data = subset(pdat,!is.na(cs)),mapping = aes(color = cs),
             shape = 1,size = 1.5) +
  geom_text_repel(color = "black",size = 2.25,min.segment.length = 0,
                  max.overlaps = Inf) +
  scale_color_manual(values = cs_colors,na.value = "darkgray") +
  guides(color = guide_legend(override.aes = list(shape = 20,size = 2))) +
  labs(x = "base-pair position on chromosome 1 (Mb)",y = "PIP") +
  theme_cowplot(font_size = 10)
print(p)
