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
# Colors are from colorbrewer2.org
cs_colors <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99",
               "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#b15928",
               "#8dd3c7","#bebada","#fb8072","#80b1d3","#fdb462",
               "#b3de69","#fccde5","#bc80bd")
zoomin_region <- c(56.34,56.55)

# Import the Illumina HumanMethylation450 CpG probe data.
load("../data/IlluminaHumanMethylation450K.rda")
probes <- IlluminaHumanMethylation450K
probes <- subset(probes,!is.na(Start_hg38) & !is.na(End_hg38))
probes <- subset(probes,End_hg38 > Start_hg38)
probes <- subset(probes,CHR_hg38 == "chr6")
probes <- probes[c("CHR_hg38","Start_hg38","End_hg38")]
names(probes) <- c("chr","start","end")
rownames(probes) <- NULL

# TOP PANEL: PIP PLOT
# -------------------
load("../outputs/CASS4_obj.RData")
ids <- names(obj_plot$fsusie_obj_me$pip)
pdat <- data.frame(id        = as.character(NA),
                   pos       = sapply(strsplit(ids,":",fixed = TRUE),"[[",2),
                   pip       = obj_plot$fsusie_obj_me$pip,
                   cs        = as.numeric(NA),
                   cs_label  = as.character(NA),
                   blacklist = FALSE,
                   stringsAsFactors = FALSE)
rownames(pdat) <- ids
n <- length(obj_plot$fsusie_obj_me$sets$cs)
for (i in 1:n) {
  snps <- names(obj_plot$fsusie_obj_me$sets$cs[[i]])
  j <- snps[which.max(pdat[snps,"pip"])]
  pdat[snps,"cs"] <- i
  pdat[snps,"cs_label"] <- sprintf("CS %d (%d SNPs, %s)",i,length(snps),j)
  pdat[j,"id"] <- paste("CS",i)
}

# Remove CSs with SNPs inside CpG probes.
snps <- which(!is.na(pdat$cs))
for (i in snps) {
  pos <- pdat[i,"pos"]
  if (any(pos >= probes$start &
          pos <= probes$end)) {
     pdat[i,"blacklist"] <- TRUE
   }
}
blacklisted_cs <- unique(subset(pdat,blacklist)$cs)
pdat <- subset(pdat,!is.element(cs,blacklisted_cs))
pdat <- transform(pdat,
                  cs       = factor(cs),
                  cs_label = factor(cs_label),
                  pos      = as.numeric(pos)/1e6)
zoomout_region <- range(pdat$pos)

# Thin out the SNPs with very small PIPs; visually, they are not all
# needed.
i <- which(pdat$pip >= 0.01)
j <- which(pdat$pip < 0.01)
j <- sample(j,2400)
rows <- sort(c(i,j))
pdat <- pdat[rows,]

# Create the PIP plot.
p1 <- ggplot(pdat,aes(x = pos,y = pip,label = id)) +
  geom_point(color = "black",size = 0.65) +
  geom_point(data = subset(pdat,!is.na(cs)),mapping = aes(color = cs_label)) +
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
  guides(color = guide_legend(nrow = 6)) +
  labs(x = "base-pair position on chromosome 20 (Mb)",y = "PIP",
       title = "CASS4",color = "") +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "top")

# BOTTOM PANEL: EFFECT PLOT
# -------------------------
source("../scripts_plot/cases_study/interpolate_effect_estimates.R")
load("../outputs/CASS4_all_effects.RData")
keep_cs <- as.numeric(levels(pdat$cs))
effects <- NULL
for (i in keep_cs) {
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
                   cs  = paste("CS",i),
                   nonzero = !(low <= 0 & up >= 0))
  j <- which.max(abs(dat$effect))
  dat[j,"label"] <- paste("CS",i)
  effects <- rbind(effects,dat)
}
effects <- transform(effects,
                     cs = factor(cs),
                     effect = effect/abs(up - low))
p2 <- ggplot() +
  geom_hline(yintercept = 0,linetype = "dotted") +
  geom_point(data = effects,
             mapping = aes(x = pos,y = effect,color = cs)) +
  geom_errorbarh(data = data.frame(xmin = zoomin_region[1],
                                   xmax = zoomin_region[2],
                                   y = -0.5),
                 mapping = aes(xmin = xmin,xmax = xmax,y = y),
                 color = "darkgray",linewidth = 0.5,height = 0.03,
                 inherit.aes = FALSE) +
  geom_text_repel(data = effects,
                  mapping = aes(x = pos,y = effect,label = label),
                  color = "black",size = 2.25,min.segment.length = 0,
                  max.overlaps = Inf,segment.color = "black") +
  scale_x_continuous(limits = zoomout_region,
                     breaks = seq(200,210,0.5)) +
  scale_color_manual(values = cs_colors,na.value = "darkgray",
                     drop = FALSE) +
  guides(color = guide_legend(nrow = 2)) +
  labs(x = "base-pair position on chromosome 20 (Mb)",y = "effect",
       color = "") +
  theme_cowplot(font_size = 8) +
  theme(legend.position = "bottom")

# Create the final combined plot.
print(plot_grid(p1,p2,nrow = 2,ncol = 1,rel_heights = c(10,6),align = "v"))
ggsave("zoomout_cass4.pdf",
       plot_grid(p1,p2,nrow = 2,ncol = 1,rel_heights = c(10,6),align = "v"),
       height = 4,width = 6)
