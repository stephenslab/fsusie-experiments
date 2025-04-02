# Create a PIP "zoomout" plot for the CASS4 locus.
#
# NOTE: download CASS4_obj.RData from
# https://uchicago.box.com/s/tt1vgg7vqayfthg0vsbl8gw0sdi0f1uo
library(data.table)
library(ggplot2)
library(ggrepel)
library(cowplot)
set.seed(1)
zoomin_region <- c(56.3,56.6)
cs_colors <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99",
               "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#b15928",
               "#8dd3c7","#bebada","#fb8072","#80b1d3","#fdb462",
               "#b3de69","#fccde5","#bc80bd")
load("../outputs/CASS4_obj.RData")
ids <- names(obj_plot$fsusie_obj_me$pip)
pdat <- data.frame(id  = as.character(NA),
                   pos = sapply(strsplit(ids,":",fixed = TRUE),"[[",2),
                   pip = obj_plot$fsusie_obj_me$pip,
                   cs  = as.character(NA),
                   stringsAsFactors = FALSE)
pdat <- transform(pdat,pos = as.numeric(pos)/1e6)
rownames(pdat) <- ids
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
j <- sample(j,2500)
i <- sort(c(i,j))
pdat <- pdat[i,]
p <- ggplot(pdat,aes(x = pos,y = pip,label = id)) +
  geom_point(color = "black",size = 0.5) +
  geom_point(data = subset(pdat,!is.na(cs)),mapping = aes(color = cs),
             shape = 1,size = 1.75) +
  geom_text_repel(color = "dimgray",size = 2.25,min.segment.length = 0,
                  max.overlaps = Inf,segment.color = "dimgray") +
  geom_errorbarh(data = data.frame(xmin = zoomin_region[1],
                                   xmax = zoomin_region[2],
                                   y = -0.1),
               mapping = aes(xmin = xmin,xmax = xmax,y = y),
               color = "darkgray",linewidth = 0.5,height = 0.1,
               inherit.aes = FALSE) +
  scale_x_continuous(breaks = seq(53,58,0.5)) +
  scale_color_manual(values = cs_colors,na.value = "darkgray") +
  ylim(-0.2,2.5) +
  # guides(color = guide_legend(override.aes = list(shape = 20,size = 2))) +
  guides(color = "none") +
  labs(x = "base-pair position on chromosome 20 (Mb)",y = "PIP",
       title = "CASS4") +
  theme_cowplot(font_size = 8)
print(p)
ggsave("zoomout_cass4.pdf",p,height = 2,width = 8)
