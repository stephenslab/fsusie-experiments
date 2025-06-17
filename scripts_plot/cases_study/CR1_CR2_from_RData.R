# Create a "zoom-in" plot for the CR1/CR2 locus.
#
# NOTE: download CR1_CR2_obj.RData from
# https://uchicago.box.com/s/tt1vgg7vqayfthg0vsbl8gw0sdi0f1uo
library(data.table)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(cowplot)
source("get_gene_annotations.R")
source("interpolate_effect_estimates.R")
load("../../outputs/CR1_CR2_obj.RData")
gene_file <-
  file.path("../../data/genome_annotations",
    "Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz")
genes <- get_gene_annotations(gene_file)

pos0 <- 207.45e6
pos1 <- 207.75e6
key_marker <- 207.577223

# The top panel shows the Alzheimer's Disease (AD) association
# p-values.
pdat1 <- as.data.frame(obj_plot$plot_df)
pdat1 <- subset(pdat1,
                study == "AD_Bellenguez_2022" &
                pos >= pos0 &
                pos <= pos1)
pdat1 <- pdat1[c("variant_alternate_id","pos","CS1","-log10(P)")]
pdat1 <- transform(pdat1,pos = pos/1e6)
names(pdat1) <- c("id","pos","CS","pval")
# subset(pdat1,id == "chr1:207577223:T:C")
#                 id   pos   CS  pval
# chr1:207577223:T:C 207.6 TRUE 32.25
ids <- pdat1$id
ids[] <- NA
ids[pdat1$id == "chr1:207510847:T:G"] <- "rs12037841"
ids[pdat1$id == "chr1:207577223:T:C"] <- "rs679515"
ids[pdat1$id == "chr1:207623552:A:T"] <- "rs10863417"
ids[pdat1$id == "chr1:207624893:C:G"] <- "rs10863418"
ids[pdat1$id == "chr1:207629207:A:C"] <- "rs4844610"
pdat1$id <- ids
# > subset(pdat1,CS)
#                 id   pos   CS  pval
# chr1:207510847:T:G 207.5 TRUE 31.75
# chr1:207577223:T:C 207.6 TRUE 32.25
# chr1:207623552:A:T 207.6 TRUE 30.59
# chr1:207624893:C:G 207.6 TRUE 30.67
# chr1:207629207:A:C 207.6 TRUE 31.95
p1 <- ggplot(pdat1,aes(x = pos,y = pval,color = CS,label = id)) +
  geom_point(size = 0.75) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  geom_text_repel(size = 2.25,color = "dimgray",segment.color = "dimgray",
                  min.segment.length = 0,max.overlaps = Inf) +
  scale_color_manual(values = c("black","dodgerblue")) +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(207,208,0.05),
                     labels = NULL) +
  ylim(0,45) +
  labs(x = "",y = "AD") + 
  theme_cowplot(font_size = 9)

# The second panel shows the CR1 eQTL p-values.
pdat2 <- as.data.frame(obj_plot$plot_df)
pdat2 <- subset(pdat2,
                study == "DLPFC_DeJager_eQTL" &
                region == "CR1" &
                pos >= pos0 &
                pos <= pos1)
pdat2$CS <- as.numeric(NA)
pdat2[which(pdat2$CS1),"CS"] <- 1
pdat2[which(pdat2$CS3),"CS"] <- 3
pdat2 <- pdat2[c("variant_alternate_id","pos","CS","-log10(P)","z")]
names(pdat2) <- c("id","pos","CS","pval","z")
pdat2 <- transform(pdat2,
                   pos = pos/1e6,
                   CS  = factor(CS))
# > subset(pdat2,id == "chr1:207577223:T:C")
#                 id   pos CS  pval      z
# chr1:207577223:T:C 207.6  1 64.89 -17.11
#
# > subset(pdat2,CS == 1)
#                 id   pos   CS  pval
# chr1:207564732:T:C 207.6 TRUE 64.56
# chr1:207573951:A:G 207.6 TRUE 63.47
# chr1:207577223:T:C 207.6 TRUE 64.89
# chr1:207611623:A:G 207.6 TRUE 64.47
# chr1:207612944:A:G 207.6 TRUE 64.47
# chr1:207613197:A:G 207.6 TRUE 64.51
# chr1:207613483:A:G 207.6 TRUE 64.47
# chr1:207629207:A:C 207.6 TRUE 62.42
# > nrow(subset(pdat2,CS == 1))
# 8
# > nrow(subset(pdat2,CS == 3))
# 12
#
# > intersect(subset(pdat1,CS)$pos,subset(pdat2,CS == 1)$pos)*1e6
# 207577223 207629207
ids <- pdat2$id
ids[] <- NA
ids[pdat2$id == "chr1:207577223:T:C"] <- "rs679515"
pdat2$id <- ids
p2 <- ggplot(pdat2,aes(x = pos,y = pval,color = CS,label = id)) +
  geom_point(size = 0.75) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  scale_color_manual(values = c("limegreen","gold"),na.value = "black") +
  geom_text_repel(size = 2.25,color = "dimgray",segment.color = "dimgray",
                  min.segment.length = 0,max.overlaps = Inf) +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(207,208,0.05),
                     labels = NULL) +
  labs(x = "",y = "CR1 eQTL in DLPFC") + 
  theme_cowplot(font_size = 9)

# The third panel shows the CR2 eQTL p-values.
pdat3 <- as.data.frame(obj_plot$plot_df)
pdat3 <- subset(pdat3,
                study == "DLPFC_DeJager_eQTL" &
                region == "CR2" &
                pos >= pos0 &
                pos <= pos1)
pdat3 <- pdat3[c("variant_alternate_id","pos","CS1","-log10(P)","z")]
names(pdat3) <- c("id","pos","CS","pval","z")
pdat3 <- transform(pdat3,pos = pos/1e6)
# > subset(pdat3,id == "chr1:207577223:T:C")
#                 id   pos   CS  pval      z
# chr1:207577223:T:C 207.6 TRUE 9.661 -6.348
# > nrow(subset(pdat3,CS))
# 20
# > intersect(subset(pdat1,CS)$pos,subset(pdat3,CS)$pos)*1e6
# 207510847 207577223
ids <- pdat3$id
ids[] <- NA
ids[pdat3$id == "chr1:207577223:T:C"] <- "rs679515"
pdat3$id <- ids
p3 <- ggplot(pdat3,aes(x = pos,y = pval,color = CS,label = id)) +
  geom_point(size = 0.75) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  scale_color_manual(values = c("black","darkorchid"),na.value = "black") +
  geom_text_repel(size = 2.25,color = "dimgray",segment.color = "dimgray",
                  min.segment.length = 0,max.overlaps = Inf) +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(207,208,0.05),
                     labels = NULL) +
  labs(x = "",y = "CR2 eQTL in DLPFC") + 
  theme_cowplot(font_size = 9)

# The fourth panel shows the haSNP PIPs.
ids   <- names(obj_plot$pip_fsusie_obj)
pdat4 <- data.frame(id  = as.character(NA),
                    pos = sapply(strsplit(ids,":",fixed = TRUE),"[[",2),
                    pip = obj_plot$pip_fsusie_obj,
                    cs  = as.character(NA),
                    stringsAsFactors = FALSE)
pdat4 <- transform(pdat4,pos = as.numeric(pos))
n <- length(obj_plot$cs_fsusie_obj)
for (i in 1:n) {
  snps <- names(obj_plot$cs_fsusie_obj[[i]])
  pdat4[snps,"cs"] <- i
}
pdat4 <- subset(pdat4,pos >= pos0 & pos <= pos1)
pdat4 <- transform(pdat4,pos = pos/1e6)
# > nrow(subset(pdat4,cs == 5))
# 15
# > subset(pdat4,cs == 5 & pip > 0.05)
#                           id   pos    pip cs
# chr1:207577223:T:C      <NA> 207.6 0.1375  5
# chr1:207598421:CT:CTT   <NA> 207.6 0.3343  5
# chr1:207619376:CAAA:CAA <NA> 207.6 0.2049  5
pdat4[c("chr1:207577223:T:C",
        "chr1:207598421:CT:CTT",
        "chr1:207619376:CAAA:CAA"),"id"] <-
  c("rs679515","rs11392366","rs869302047")
p4 <- ggplot(pdat4,aes(x = pos,y = pip,color = cs,label = id)) +
  geom_point(size = 0.75) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  geom_text_repel(size = 2.25,color = "dimgray",segment.color = "dimgray",
                  min.segment.length = 0,max.overlaps = Inf) +
  scale_color_manual(values = c("tomato","cyan"),na.value = "black") +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(207,208,0.05),
                     labels = NULL) +
  ylim(0,0.4) +
  labs(x = "",y = "haSNP PIP") + 
  theme_cowplot(font_size = 9)

# The fifth panel shows the estimated effects on the H3K27ac levels.
pdat5 <- data.frame(effect = obj_plot$effect_s[1,],
                    up     = obj_plot$effect_s[2,],
                    low    = obj_plot$effect_s[3,],
                    pos    = obj_plot$pos_H3Kac_effect)
n     <- length(obj_plot$peak_pos)
pdat5 <- interpolate_effect_estimates(pdat5,obj_plot$peak_pos[seq(2,n-1)])
pdat5 <- subset(pdat5,
                pos >= pos0 &
                pos <= pos1)
pdat5 <- transform(pdat5,
                   pos = pos/1e6,
                   y   = 0)
pdat5 <- list(zero_effects    = subset(pdat5,low <= 0 & up >= 0),
              nonzero_effects = subset(pdat5,!(low <= 0 & up >= 0)))
p5 <- ggplot() +
  geom_hline(yintercept = 0,linetype = "dotted") +
  geom_linerange(data = pdat5$nonzero_effects,
                 mapping = aes(x = pos,ymin = low,ymax = up),
                 color = "dodgerblue") +
  geom_point(data = pdat5$nonzero_effects,
             mapping = aes(x = pos,y = effect),
             shape = 21,color = "white",fill = "dodgerblue",size = 1.1) +
  geom_point(data = pdat5$zero_effects,
             mapping = aes(x = pos,y = y),
             shape = 21,color = "white",fill = "dodgerblue",size = 1.1) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(207,208,0.05),
                     labels = NULL) +
  labs(x = "",y = "effect of C allele") +
  theme_cowplot(font_size = 9)

# The sixth panel shows the raw data.
pdat6 <- obj_plot$count_df
names(pdat6) <- c("CC","CT","TT","pos")
pdat6 <- subset(pdat6,pos >= pos0 & pos <= pos1)
rows1 <- with(pdat6,which(pmax(TT,CT,CC) >= 5))
rows2 <- with(pdat6,which(pmax(TT,CT,CC) < 5))
rows2 <- sample(rows2,5000)
rows  <- c(rows1,rows2)
pdat6 <- pdat6[rows,]
pdat6 <- transform(pdat6,pos = pos/1e6)
pdat6 <- reshape2::melt(pdat6,id.vars = "pos",variable.name = "genotype",
                        value.name = "count")
pdat6 <- transform(pdat6,genotype = factor(genotype))
rows  <- order(pdat6$genotype,decreasing = TRUE)
pdat6 <- pdat6[rows,]
p6 <- ggplot(pdat6,aes(x = pos,y = count,color = genotype)) +
  geom_point(size = 0.35) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  scale_color_manual(values = c("darkblue","darkviolet","darkorange")) +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(207,208,0.05),
                     labels = NULL) +
  ylim(2,24) +
  labs(x = "") + 
  theme_cowplot(font_size = 9)
             
# The seventh panel shows the genes.
pdat7 <- subset(genes,
                chromosome == "chr1" &
                end > pos0 &
                start < pos1)
pdat7 <- transform(pdat7,tss = ifelse(strand == "+",start,end))
pdat7 <- transform(pdat7,
                   start = start/1e6,
                   end   = end/1e6,
                   tss   = tss/1e6)
n <- nrow(pdat7)
pdat7$y <- seq(0,1,length.out = n)
p7 <- ggplot(pdat7,aes(x = start,xend = end,y = y,yend = y,
                       label = gene_name)) +
  geom_segment(color = "dodgerblue",linewidth = 0.5) +
  geom_point(mapping = aes(x = tss),color = "dodgerblue",size = 1.5,
             shape = 18) +
  geom_text(color = "black",size = 2.25,fontface = "italic",
            hjust = "right",nudge_x = -0.003) +
  geom_vline(xintercept = key_marker,linetype = "dotted",color = "darkgray") +
  scale_x_continuous(limits = c(pos0,pos1)/1e6,
                     breaks = seq(207,208,0.05)) +
  scale_y_continuous(limits = c(-0.1,1.1),breaks = NULL) +
  labs(x = "base-pair position on chromosome 1 (Mb)",y = "") + 
  theme_cowplot(font_size = 9)

# Save the full figure to a PDF.
print(plot_grid(p1,p2,p3,p4,p5,p6,p7,nrow = 7,ncol = 1,align = "v"))
ggsave("CR1CR2_zoomin_plot.pdf",
       plot_grid(p1,p2,p3,p4,p5,p6,p7,nrow = 7,ncol = 1,align = "v"),
       height = 5.75,width = 4.5)
