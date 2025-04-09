# Create a "zoom-in" plot for the CD2AP locus.
#
# NOTE: download CD2AP_obj.RData from
# https://uchicago.box.com/s/tt1vgg7vqayfthg0vsbl8gw0sdi0f1uo
library(data.table)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(cowplot)
source("get_gene_annotations.R")
source("interpolate_effect_estimates.R")
gene_file <-
  file.path("../../data/genome_annotations",
    "Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz")
genes <- get_gene_annotations(gene_file)
load("../../outputs/CD2AP_obj.RData")

pos0 <- 47.375e6
pos1 <- 47.725e6
key_marker <- 47472829/1e6

# The top panel shows the Alzheimer's Disease (AD) association
# p-values.
pdat1 <- as.data.frame(obj_plot$pdat)
pdat1 <- subset(pdat1,
                study == "AD_Bellenguez_2022" &
                pos >= pos0 &
                pos <= pos1)
pdat1 <- pdat1[c("variant_alternate_id","pos","CS1","-log10(P)","z")]
pdat1 <- transform(pdat1,pos = pos/1e6)
names(pdat1) <- c("id","pos","CS","pval","z")
AD_snps <- subset(pdat1,CS)$id
# > subset(pdat1,id == "chr6:47472829:C:A")
#                id   pos   CS  pval     z
# chr6:47472829:C:A 47.47 TRUE 11.22 6.879
# > nrow(subset(pdat1,CS))
# 34
ids <- pdat1$id
ids[] <- NA
ids[pdat1$id == "chr6:47472829:C:A"] <- "rs9369695"
ids[pdat1$id == "chr6:47615928:C:T"] <- "rs12195738"
pdat1$id <- ids
p1 <- ggplot(pdat1,aes(x = pos,y = pval,color = CS,label = id)) +
  geom_point(size = 0.75) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  geom_text_repel(size = 2.25,color = "dimgray",segment.color = "dimgray",
                  min.segment.length = 0,max.overlaps = Inf) +
  scale_color_manual(values = c("black","dodgerblue")) +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(45,48,0.05),
                     labels = NULL) +
  ylim(0,15) +
  labs(x = "",y = "AD") + 
  theme_cowplot(font_size = 9)

# The second panel shows the CD2AP eQTL p-values.
pdat2 <- as.data.frame(obj_plot$pdat)
pdat2 <- subset(pdat2,
                region == "CD2AP" &
                study == "Mic_DeJager_eQTL" &
                pos >= pos0 &
                pos <= pos1)
# > subset(pdat2,variant_alternate_id == "chr6:47472829:C:A")$z
# [1] 6.497
pdat2$CS <- as.numeric(NA)
pdat2[which(pdat2$CS1),"CS"] <- 1
pdat2 <- pdat2[c("variant_alternate_id","pos","CS","-log10(P)")]
names(pdat2) <- c("id","pos","CS","pval")
pdat2 <- transform(pdat2,
                   pos = pos/1e6,
                   CS  = factor(CS))
CD2AP_snps <- subset(pdat2,!is.na(CS))$id
# > nrow(subset(pdat2,CS == 1))
# 80
# > length(intersect(AD_snps,CD2AP_snps))
# 32
ids <- pdat2$id
ids[] <- NA
ids[pdat2$id == "chr6:47472829:C:A"] <- "rs9369695"
pdat2$id <- ids
p2 <- ggplot(pdat2,aes(x = pos,y = pval,color = CS,label = id)) +
  geom_point(size = 0.75) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  scale_color_manual(values = c("limegreen","gold"),na.value = "black") +
  geom_text_repel(size = 2.25,color = "dimgray",segment.color = "dimgray",
                  min.segment.length = 0,max.overlaps = Inf) +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(45,48,0.05),
                     labels = NULL) +
  ylim(0,15) +
  labs(x = "",y = "CD2AP eQTL in microglia") + 
  theme_cowplot(font_size = 9)

# The third panel shows the mSNP PIPs.
ids   <- names(obj_plot$fsusie_obj_me$pip)
pdat3 <- data.frame(id  = as.character(NA),
                    pos = sapply(strsplit(ids,":",fixed = TRUE),"[[",2),
                    pip = obj_plot$fsusie_obj_me$pip,
                    cs  = as.character(NA),
                    stringsAsFactors = FALSE)
pdat3 <- transform(pdat3,pos = as.numeric(pos))
n <- length(obj_plot$fsusie_obj_me$sets$cs)
for (i in 1:n) {
  snps <- names(obj_plot$fsusie_obj_me$sets$cs[[i]])
  pdat3[snps,"cs"] <- i
}
pdat3 <- subset(pdat3,pos >= pos0 & pos <= pos1)
pdat3 <- transform(pdat3,pos = pos/1e6)
# > subset(pdat3,!is.na(cs))
#                     id   pos    pip cs
# chr6:47405566:G:A <NA> 47.41 0.5000 16
# chr6:47405711:A:G <NA> 47.41 0.5000 16
# chr6:47472829:C:A <NA> 47.47 0.9645 13
# chr6:47664375:T:C <NA> 47.66 0.3333  8
# chr6:47664664:A:T <NA> 47.66 0.3333  8
# chr6:47665074:C:T <NA> 47.67 0.3333  8
pdat3[c("chr6:47405566:G:A",
        "chr6:47472829:C:A",
        "chr6:47664375:T:C"),"id"] <-
    c("rs4715009","rs9369695","rs5013495")
p3 <- ggplot(pdat3,aes(x = pos,y = pip,color = cs,label = id)) +
  geom_point(size = 0.75) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  geom_text_repel(size = 2.25,color = "dimgray",segment.color = "dimgray",
                  min.segment.length = 0,max.overlaps = Inf) +
  scale_color_manual(values = c("magenta","darkorange","darkviolet"),
                     na.value = "black") +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(45,48,0.05),
                     labels = NULL) +
  labs(x = "",y = "mSNP PIP") + 
  theme_cowplot(font_size = 9)

# The fourth panel shows the estimated effects on the methylation levels.
pdat4 <- data.frame(effect = obj_plot$effect_s[1,],
                    up     = obj_plot$effect_s[2,],
                    low    = obj_plot$effect_s[3,],
                    pos    = obj_plot$pos_est_effect)
n     <- length(obj_plot$me_pos)
pdat4 <- interpolate_effect_estimates(pdat4,obj_plot$me_pos[seq(2,n-1)])
pdat4 <- subset(pdat4,
                pos >= pos0 &
                pos <= pos1)
pdat4 <- transform(pdat4,
                   pos = pos/1e6,
                   y   = 0)
pdat4 <- list(zero_effects    = subset(pdat4,low <= 0 & up >= 0),
              nonzero_effects = subset(pdat4,!(low <= 0 & up >= 0)))
p4 <- ggplot() +
  geom_hline(yintercept = 0,linetype = "dotted") +
  geom_linerange(data = pdat4$nonzero_effects,
                 mapping = aes(x = pos,ymin = low,ymax = up),
                 color = "dodgerblue") +
  geom_point(data = pdat4$nonzero_effects,
             mapping = aes(x = pos,y = effect),
             color = "dodgerblue",size = 0.75) +
  geom_point(data = pdat4$zero_effects,
             mapping = aes(x = pos,y = y),
             color = "dodgerblue",size = 0.75) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(45,48,0.05),
                     labels = NULL) +
  labs(x = "",y = "effect of A allele") +
  theme_cowplot(font_size = 9)

# The fifth panel shows the genes.
pdat5 <- subset(genes,
                chromosome == "chr6" &
                end > pos0 &
                start < pos1)
pdat5 <- transform(pdat5,tss = ifelse(strand == "+",start,end))
pdat5 <- transform(pdat5,
                   start = start/1e6,
                   end   = end/1e6,
                   tss   = tss/1e6)
n <- nrow(pdat5)
pdat5$y <- seq(0,1,length.out = n)
p5 <- ggplot(pdat5,aes(x = start,xend = end,y = y,yend = y,
                       label = gene_name)) +
  geom_segment(color = "dodgerblue",linewidth = 0.5) +
  geom_point(mapping = aes(x = tss),color = "dodgerblue",size = 1.5,
             shape = 18) +
  geom_text(color = "black",size = 2.25,fontface = "italic",
            hjust = "right",nudge_x = -0.003) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(45,48,0.05)) +
  scale_y_continuous(limits = c(-0.1,1.1),breaks = NULL) +
  labs(x = "base-pair position on chromosome 6 (Mb)",y = "") + 
  theme_cowplot(font_size = 9)

# Save the full figure to a PDF.
print(plot_grid(p1,p2,p3,p4,p5,nrow = 5,ncol = 1,align = "v"))
ggsave("CD2AP_zoomin_plot.pdf",
       plot_grid(p1,p2,p3,p4,p5,nrow = 5,ncol = 1,align = "v"),
       height = 4.25,width = 4.5)

