# Create a "zoom-in" plot for the CASS4 locus.
#
# NOTE: download CASS4_obj.RData from
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
load("../../outputs/CASS4_obj.RData")

pos0 <- 56.34e6
pos1 <- 56.55e6
key_marker <- 56409008/1e6

# The top panel shows the Alzheimer's Disease (AD) association
# p-values.
lpfromz <- function (z)
  (log(2) + pnorm(-abs(z),log.p = TRUE))/log(10)
pdat1 <- as.data.frame(obj_plot$AD_GWAS)
pdat1 <- transform(pdat1,
                   pval = -lpfromz(z))
pdat1 <- subset(pdat1,
                pos >= pos0 &
                pos <= pos1)
pdat1 <- pdat1[c("variant_alternate_id","pos","CS1","pval","z")]
pdat1 <- transform(pdat1,pos = pos/1e6)
names(pdat1) <- c("id","pos","CS","pval","z")
AD_snps <- subset(pdat1,CS)$id
# > subset(pdat1,id == "chr20:56409008:G:C")
#                 id   pos   CS  pval
# chr20:56409008:G:C 56.41 TRUE 8.053
# > nrow(subset(pdat1,CS))
# 13
ids <- pdat1$id
ids[] <- NA
ids[pdat1$id == "chr20:56409008:G:C"] <- "rs1884913"
ids[pdat1$id == "chr20:56423488:A:G"] <- "rs6014724"
pdat1$id <- ids
p1 <- ggplot(pdat1,aes(x = pos,y = pval,color = CS,label = id)) +
  geom_point(size = 0.75) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  geom_text_repel(size = 2.25,color = "dimgray",segment.color = "dimgray",
                  min.segment.length = 0,max.overlaps = Inf) +
  scale_color_manual(values = c("black","dodgerblue")) +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(56,57,0.05),
                     labels = NULL) +
  ylim(0,12) +
  labs(x = "",y = "AD") + 
  theme_cowplot(font_size = 9)

# The second panel shows the CASS4 eQTL p-values.
pdat2 <- as.data.frame(obj_plot$qTLdata)
pdat2 <- transform(pdat2,
                   pval = -lpfromz(z))
pdat2 <- subset(pdat2,
                region == "CASS4" &
                study == "DLPFC_DeJager_eQTL" &
                pos >= pos0 &
                pos <= pos1)
pdat2 <- pdat2[c("variant_id","pos","cs_coverage_0.95","pval")]
names(pdat2) <- c("id","pos","CS","pval")
pdat2 <- transform(pdat2,
                   pos = pos/1e6,
                   CS  = factor(CS))
CASS4_snps <- subset(pdat2,!is.na(CS))$id
# > nrow(subset(pdat2,CS == 0))
# 11
# > length(intersect(AD_snps,CASS4_snps))
# 0
ids <- pdat2$id
ids[] <- NA
ids[pdat2$id == "chr20:56438160:A:G"] <- "rs6014730"
ids[pdat2$id == "chr20:56409008:G:C"] <- "rs1884913"
pdat2$id <- ids
p2 <- ggplot(pdat2,aes(x = pos,y = pval,color = CS,label = id)) +
  geom_point(size = 0.75) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  scale_color_manual(values = c("limegreen","gold"),na.value = "black") +
  geom_text_repel(size = 2.25,color = "dimgray",segment.color = "dimgray",
                  min.segment.length = 0,max.overlaps = Inf) +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(56,57,0.05),
                     labels = NULL) +
  labs(x = "",y = "CASS4 eQTL in DLPFC") + 
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
pdat3 <- transform(pdat3,
                   pos = pos/1e6,
                   cs = factor(cs))
# > subset(pdat3,cs == 4)
#                      id   pos    pip cs
# chr20:56348142:G:A <NA> 56.35 0.3318  4
# chr20:56348706:T:C <NA> 56.35 0.3318  4
# chr20:56348928:C:T <NA> 56.35 0.3318  4
# > subset(pdat3,cs == 15)
#                        id   pos    pip cs
# chr20:56393188:G:GTA <NA> 56.39 0.8602 15
# chr20:56399584:A:G   <NA> 56.40 0.1131 15
# > subset(pdat3,cs == 16)
#                      id   pos pip cs
# chr20:56379506:T:C <NA> 56.38   1 16
# > subset(pdat3,cs == 18)
#                      id   pos    pip cs
# chr20:56409008:G:C <NA> 56.41 0.9999 18
pdat3[c("chr20:56348142:G:A",
        "chr20:56393188:G:GTA",
        "chr20:56379506:T:C",
        "chr20:56409008:G:C"),"id"] <-
    c("rs8118732","chr20:56393188","rs112535545","rs1884913")
p3 <- ggplot(pdat3,aes(x = pos,y = pip,color = cs,label = id)) +
  geom_point(size = 0.75) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  geom_text_repel(size = 2.25,color = "dimgray",segment.color = "dimgray",
                  min.segment.length = 0,max.overlaps = Inf) +
  scale_color_manual(values = c("magenta","darkorange","forestgreen",
                                "darkviolet"),
                     na.value = "black") +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(56,57,0.05),
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
                     breaks = seq(56,57,0.05),
                     labels = NULL) +
  labs(x = "",y = "effect of C allele") +
  theme_cowplot(font_size = 9)

# The fifth panel shows the genes.
pdat5 <- subset(genes,
                chromosome == "chr20" &
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
                     breaks = seq(56,57,0.05)) +
  scale_y_continuous(limits = c(-0.1,1.1),breaks = NULL) +
  labs(x = "base-pair position on chromosome 20 (Mb)",y = "") + 
  theme_cowplot(font_size = 9)

# Save the full figure to a PDF.
print(plot_grid(p1,p2,p3,p4,p5,nrow = 5,ncol = 1,align = "v"))
ggsave("CASS4_zoomin_plot.pdf",
       plot_grid(p1,p2,p3,p4,p5,nrow = 5,ncol = 1,align = "v"),
       height = 4.25,width = 4.5)
