# Create a PIP "zoomout" plot for the CR1/CR2 locus.
library(data.table)
library(ggplot2)
library(ggrepel)
library(cowplot)
zoomin_region <- c(207.4,207.7)
# Colors are from colorbrewer2.org
cs_colors <- c("dodgerblue","darkorange","magenta","limegreen","gold")
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
  cs_label <- sprintf("(CS %d, %d SNPs)",i,length(snps))
  pdat[snps,"cs"] <- i
  j <- snps[which.max(pdat[snps,"pip"])]
  pdat[j,"id"] <- paste(j,cs_label)
}
pdat <- transform(pdat,cs = factor(cs))
p <- ggplot(pdat,aes(x = pos,y = pip,label = id)) +
  geom_point(color = "black",size = 0.5) +
  geom_point(data = subset(pdat,!is.na(cs)),mapping = aes(color = cs),
             shape = 1,size = 1.5) +
  geom_text_repel(color = "black",size = 2.25,min.segment.length = 0,
                  max.overlaps = Inf) +
  geom_errorbarh(data = data.frame(xmin = zoomin_region[1],
                                   xmax = zoomin_region[2],
                                   y = -0.1),
                 mapping = aes(xmin = xmin,xmax = xmax,y = y),
                 color = "darkgray",linewidth = 0.5,height = 0.05,
                 inherit.aes = FALSE) +
  scale_color_manual(values = cs_colors,na.value = "darkgray") +
  guides(color = "none") +
# guides(guide_legend(override.aes = list(shape = 20,size = 2))) +
  labs(x = "base-pair position on chromosome 1 (Mb)",y = "PIP",
       title = "CR1/CR2") +
  theme_cowplot(font_size = 10)
print(p)
ggsave("zoomout_cr1.pdf",p,height = 2,width = 8)
