# Create a PIP "zoomout" plot for the CD2AP locus.
#
# NOTE: download CD2AP_obj.RData and CD2AP_all_effects.RData from
# https://uchicago.box.com/s/tt1vgg7vqayfthg0vsbl8gw0sdi0f1uo
#
# SEE ALSO:
# scripts_plot/cases_study/CD2AP_from_RData.R
#
library(data.table)
library(ggplot2)
library(ggrepel)
library(cowplot)
set.seed(1)
zoomin_region <- c(47.375, 47.725)

# Import the Illumina HumanMethylation450 CpG probe information.
load("../data/IlluminaHumanMethylation450K.rda")
probes <- IlluminaHumanMethylation450K
probes <- subset(probes,!is.na(Start_hg38) & !is.na(End_hg38))
probes <- subset(probes,End_hg38 > Start_hg38)
probes <- probes[c("CHR_hg38","Start_hg38","End_hg38")]
names(probes) <- c("chr","start","end")
probes$cpg <- rownames(probes)
rownames(probes) <- NULL
probes <- subset(probes,chr == "chr6")

# The top panel is a PIP plot.
# Colors are from colorbrewer2.org
cs_colors <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99",
               "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#b15928",
               "#8dd3c7","#bebada","#fb8072","#80b1d3","#fdb462",
               "#b3de69","#fccde5","#bc80bd","#ccebc5","#ffed6f")
load("../outputs/CD2AP_obj.RData")
ids <- names(obj_plot$fsusie_obj_me$pip)
pdat <- data.frame(id  = as.character(NA),
                   pos = sapply(strsplit(ids,":",fixed = TRUE),"[[",2),
                   pip = obj_plot$fsusie_obj_me$pip,
                   cs  = as.character(NA),
                   cs_idx = as.numeric(NA),
                   blacklist = FALSE,
                   stringsAsFactors = FALSE)
rownames(pdat) <- ids
n <- length(obj_plot$fsusie_obj_me$sets$cs)
for (i in 1:n) {
  snps <- names(obj_plot$fsusie_obj_me$sets$cs[[i]])
  j <- snps[which.max(pdat[snps,"pip"])]
  pdat[snps,"cs_idx"] <- i
  pdat[snps,"cs"] <- sprintf("CS %d (%d SNPs, %s)",i,length(snps),j)
  pdat[j,"id"] <- paste("CS",i)
}
snps <- which(!is.na(pdat$cs))
for (i in snps) {
  pos <- pdat[i,"pos"]
  if (any(probes$start <= pos & pos <= probes$end)) {
     pdat[i,"blacklist"] <- TRUE
   }
}
blacklisted_cs <- unique(subset(pdat,blacklist)$cs_idx)
pdat <- subset(pdat,!is.element(cs_idx,blacklisted_cs))
pdat <- transform(pdat,
                  cs  = factor(cs),
                  pos = as.numeric(pos)/1e6)
zoomout_region <- range(pdat$pos)
i <- which(pdat$pip >= 0.01)
j <- which(pdat$pip < 0.01)
j <- sample(j,2000)
i <- sort(c(i,j))
pdat <- pdat[i,]
p1 <- ggplot(pdat,aes(x = pos,y = pip,label = id)) +
  geom_point(color = "black",size = 0.65) +
  geom_point(data = subset(pdat,!is.na(cs)),mapping = aes(color = cs)) +
  geom_text_repel(color = "black",size = 2.25,min.segment.length = 0,
                  max.overlaps = Inf,segment.color = "black") +
  geom_errorbarh(data = data.frame(xmin = zoomin_region[1],
                                   xmax = zoomin_region[2],
                                   y = -0.1),
                 mapping = aes(xmin = xmin,xmax = xmax,y = y),
                 color = "darkgray",linewidth = 0.5,height = 0.05,
                 inherit.aes = FALSE) +
  scale_x_continuous(breaks = seq(44,49,0.5)) +
  scale_color_manual(values = cs_colors,na.value = "darkgray") +
  labs(x = "base-pair position on chromosome 6 (Mb)",y = "PIP",
       title = "CD2AP") +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "top",
        legend.text = element_text(size = 6))

# Now create the plot showing the effects on methylation.
source("../scripts_plot/cases_study/interpolate_effect_estimates.R")
load("../outputs/CD2AP_all_effects.RData")
n <- length(obj_plot$effect_list)
good_cs <- setdiff(1:n,blacklisted_cs)
effects <- NULL
for (i in good_cs) {
  dat <- obj_plot$effect_list[[i]]
  me_pos <- dat$me_pos
  dat <- data.frame(effect = dat$effect_s$effect_estimate,
                    up     = dat$effect_s$cred_band["up",],
                    low    = dat$effect_s$cred_band["low",],
                    pos    = dat$pos_est_effect,
                    label  = "",
                    stringsAsFactors = FALSE)
  m   <- length(me_pos)
  dat <- interpolate_effect_estimates(dat,me_pos[seq(2,m-1)])
  dat <- transform(dat,
                   pos = pos/1e6,
                   y   = 0,
                   cs  = paste("CS",i),
                   nonzero = !(low <= 0 & up >= 0))
  j <- which.max(abs(dat$effect))
  dat[j,"label"] <- paste("CS",i)
  effects <- rbind(effects,dat)
}
effects <- transform(effects,cs = factor(cs))
effects <- transform(effects,effect = effect/abs(up - low))
p2 <- ggplot() +
  geom_hline(yintercept = 0,linetype = "dotted") +
  # geom_linerange(data = subset(effects,nonzero),
  #                mapping = aes(x = pos,ymin = low,ymax = up,color = cs)) +
  geom_point(data = effects,
             mapping = aes(x = pos,y = effect,color = cs,shape = nonzero)) +
  geom_errorbarh(data = data.frame(xmin = zoomin_region[1],
                                   xmax = zoomin_region[2],
                                   y = -20),
                 mapping = aes(xmin = xmin,xmax = xmax,y = y),
                 color = "darkgray",linewidth = 0.5,height = 1,
                 inherit.aes = FALSE) +
  geom_text_repel(data = effects,
                  mapping = aes(x = pos,y = effect,label = label),
                  color = "black",size = 2.25,min.segment.length = 0,
                  max.overlaps = Inf,segment.color = "black") +
  scale_x_continuous(limits = zoomout_region,
                     breaks = seq(44,49,0.5)) +
  scale_color_manual(values = cs_colors,na.value = "darkgray") +
  scale_shape_manual(values = c(19,19)) +
  guides(shape = "none") +
  labs(x = "base-pair position on chromosome 6 (Mb)",y = "effect") +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "top",
        legend.text = element_text(size = 6))

# Create the final combined plot.
print(plot_grid(p1,p2,nrow = 2,ncol = 1,align = "v"))

# ggsave("zoomout_cd2ap.pdf",
#        plot_grid(p1,p2,nrow = 2,ncol = 1,
#                  rel_heights = c(2,1),align = "v"),
#        height = 3,width = 8)
