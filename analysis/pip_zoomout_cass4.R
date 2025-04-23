# Create a PIP "zoomout" plot for the CASS4 locus.
#
# NOTE: download CASS4_obj.RData and CASS4_all_effects.RData from
# https://uchicago.box.com/s/tt1vgg7vqayfthg0vsbl8gw0sdi0f1uo
#
# SEE ALSO:
# scripts_plot/cases_study/CASS4_from_RData.R
#
library(data.table)
library(ggplot2)
library(ggrepel)
library(cowplot)
set.seed(1)
zoomin_region <- c(56.34,56.55)

# The top panel is a PIP plot.
# Colors are from colorbrewer2.org
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
j <- sample(j,2500)
i <- sort(c(i,j))
pdat <- pdat[i,]
p1 <- ggplot(pdat,aes(x = pos,y = pip,label = id)) +
  geom_point(color = "black",size = 0.5) +
  geom_point(data = subset(pdat,!is.na(cs)),mapping = aes(color = cs),
             shape = 1,size = 1.75) +
  geom_text_repel(color = "midnightblue",size = 2.25,min.segment.length = 0,
                  max.overlaps = Inf,segment.color = "midnightblue") +
  geom_errorbarh(data = data.frame(xmin = zoomin_region[1],
                                   xmax = zoomin_region[2],
                                   y = -0.1),
               mapping = aes(xmin = xmin,xmax = xmax,y = y),
               color = "darkgray",linewidth = 0.5,height = 0.1,
               inherit.aes = FALSE) +
  scale_x_continuous(breaks = seq(53,58,0.5)) +
  scale_color_manual(values = cs_colors,na.value = "darkgray") +
  ylim(-0.2,1.5) +
  guides(color = "none") +
  labs(x = "base-pair position on chromosome 20 (Mb)",y = "PIP",
       title = "CASS4") +
  theme_cowplot(font_size = 8)

stop()

# The second panel shows the effects on methylation.
#
# TO DO: Need to work on this code.
#
source("../scripts_plot/cases_study/interpolate_effect_estimates.R")
load("../outputs/CASS4_all_effects.RData")
n <- length(obj_plot$effect_list)
effects <- NULL
for (i in 1:n) {
  dat <- obj_plot$effect_list[[i]]
  peak_pos <- dat$peak_pos
  dat <- data.frame(effect = dat$effect_s$effect_estimate,
                    up     = dat$effect_s$cred_band["up",],
                    low    = dat$effect_s$cred_band["low",],
                    pos    = dat$pos_H3Kac_effect,
                    label  = "",
                    stringsAsFactors = FALSE)
  m   <- length(peak_pos)
  dat <- interpolate_effect_estimates(dat,peak_pos[seq(2,m-1)])
  dat <- transform(dat,
                   pos     = pos/1e6,
                   cs      = i,
                   nonzero = !(low <= 0 & up >= 0))
  j <- which.max(abs(dat$effect))
  dat[j,"label"] <- paste("CS",i)
  effects <- rbind(effects,dat)
}
effects <- transform(effects,cs = factor(cs,1:n))
p2 <- ggplot() +
  geom_hline(yintercept = 0,linetype = "dotted") +
  geom_linerange(data = subset(effects,nonzero),
                 mapping = aes(x = pos,ymin = low,ymax = up,color = cs)) +
  geom_point(data = effects,
             mapping = aes(x = pos,y = effect,color = cs,shape = nonzero),
             size = 1) +
  geom_errorbarh(data = data.frame(xmin = zoomin_region[1],
                                   xmax = zoomin_region[2],
                                   y = -0.5),
                 mapping = aes(xmin = xmin,xmax = xmax,y = y),
                 color = "darkgray",linewidth = 0.5,height = 0.03,
                 inherit.aes = FALSE) +
  geom_text_repel(data = effects,
                  mapping = aes(x = pos,y = effect,label = label),
                  color = "midnightblue",size = 2.25,min.segment.length = 0,
                  max.overlaps = Inf,segment.color = "midnightblue") +
  scale_x_continuous(limits = zoomout_region,
                     breaks = seq(200,210,0.5)) +
  scale_y_continuous(breaks = seq(-1,1,0.1)) +
  scale_color_manual(values = cs_colors,na.value = "darkgray",
                     drop = FALSE) +
  scale_shape_manual(values = c(4,19)) +
  guides(color = "none",shape = "none") +
  labs(x = "base-pair position on chromosome 20 (Mb)",y = "effect") +
  theme_cowplot(font_size = 8)

# Create the final combined plot.
# TO DO.
