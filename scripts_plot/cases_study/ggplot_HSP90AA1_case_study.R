


# HSP90aa1



rm(list=ls())

cs_colors <- c("#a6cee3","#33a02c","#1f78b4","#b2df8a","#fb9a99",
               "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#b15928",
               "#8dd3c7","#bebada","#fb8072","#80b1d3","#fdb462",
               "#b3de69","#fccde5","#bc80bd")

path="C:/Document/Serieux/Travail/Data_analysis_and_papers/GTEX_analysis_Fsusie"
source(paste0(path,"/code/plot_all_effect_log.R"))
source(paste0(path,"/code/plot_log.R"))

load(paste0(path,"/output/local/ENSG00000080824.csv.gz.RData")) 
# ---- Setup ----
library(ggplot2)
library(cowplot)
library(dplyr)


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






pdat <- transform(pdat,cs = factor(cs))
i <- which(pdat$pip >= 0.01)
j <- which(pdat$pip < 0.01)

i <- sort(c(i,j))
pdat <- pdat[i,]

min(pdat$pos)
max(pdat$pos)
p1 <- ggplot(pdat, aes(x = pos, y = pip)) +
  geom_rect(data = data.frame(xmin = zoomin_region[1],
                              xmax = zoomin_region[2],
                              ymin = -0.05,
                              ymax = 1.1),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightblue", alpha = 0.4, inherit.aes = FALSE) +
  geom_point(color = "black", size = 0.5) +
  geom_point(data = subset(pdat, !is.na(cs)), aes(color = cs),
             shape = 1, size = 1.75) +
  #geom_text_repel(
  #  data = subset(pdat, !is.na(id) & id != ""),  # Ensure valid labels
  #  mapping = aes(x = pos, y = pip, label = id),
  #  color = "midnightblue", size = 2.25,
  #  min.segment.length = 0, max.overlaps = Inf,
  #  segment.color = "midnightblue"
  #)
  
  scale_x_continuous(breaks = seq(101.8,102.2, 0.1)) +
  scale_y_continuous(limits = c(-0.05, 1.2), breaks = c(0, 0.5, 1)) +
  scale_color_manual(values = cs_colors, na.value = "darkgray") +
  guides(color = "none") +
  labs( x="", y = "PIP", title = "HSP90AA1") +
  theme_cowplot( ) +
  theme(
    axis.line.x = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(size = 9),
    #axis.title.x = element_text(size = 10, face = "bold"),
    #axis.title.y = element_text(size = 10,   face = "bold",  angle = 90, vjust = 0.5),
    #axis.text.y = element_text( vjust = 0.5, hjust = 0.5, face = "bold", size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )


p1

out$locus[1] -min(as.numeric(out$info_SNP$POS))
out$locus[2] -max(as.numeric(out$info_SNP$POS))

-200000
id1= which(  as.numeric(out$info_SNP$POS) < out$locus[1]-150000)  

library(ggtranscript)
library(ggplot2)
library(dplyr)
library(ggtranscript)
library(ggplot2)
library(dplyr)

obj= out$res
Y=log1p(as.matrix(out$Y/out$size_factor_local) )
X=as.matrix(out$X) 


obj= out$res
Y=log1p(as.matrix(out$Y/out$size_factor_local) )
X=as.matrix(out$X) 





m1 <-susieR::susie(X=X,
                   y=rowSums(Y),
                   L=5
)
m1$sets


PCA <- svd(Y)
m2 <-susieR::susie(X=X,
                   y=PCA$u[,1],
                   L=5
)

m2$sets





rest =smash_regression(obj, Y=  fsusieR:::colScale(Y , scale = FALSE)  ,
                       X= fsusieR:::colScale(X   ))

#rest= TI_regression(obj, Y=  fsusieR:::colScale(Y , scale = FALSE)  ,
#                    X= fsusieR:::colScale(X   ))
out1=out
obj1=out$res
obj=out$res
chr= paste0("chr",out1$chr)
pos0 =out1$locus[1]
pos1=out1$locus[2]

snp_info=out1$info_SNP
cs = 2
log1p_count=TRUE
data_splice=NULL
plot_cred_band=TRUE
type_data="p"

# Extract the relevant genes and exons in the specified region
region_genes <- genes(txdb,columns = c("tx_id","gene_id"))

# Subset the genes and exons to the region of interest.
region_genes <- subsetByOverlaps(region_genes,
                                 GRanges(seqnames = chr,
                                         ranges = IRanges(pos0,pos1)))

# Generate a sequence of positions with a length of 1,024.
positions <- seq(pos0,pos1,length.out = 1024)

markers <- obj$cs[[cs]]
j       <- which.max(obj$pip[markers])
marker  <- markers[j]
x       <- X[,marker]





read_counts <- rbind(colMeans(log1p(out$Y[x == 0,])),
                     colMeans(log1p(out$Y[x == 1,])) )

### CS1 -----

uni_res= univariate_functional_regression(Y= Y  ,
                                          X=as.matrix(X[,obj$cs[[cs]][1]],
                                                      ncol=1),
                                          method="TI"
)




t_effect=  smashr::smash(read_counts[2,]-read_counts[1,], post.var=TRUE,
                         sigma=sqrt(var(read_counts[2,]-read_counts[1,])))


var(read_counts[2,]-read_counts[1,])

plot(t_effect$mu.est)
cred_band = rbind( t_effect$mu.est+1.98*sqrt(t_effect$mu.est.var),
                   t_effect$mu.est-1.98*sqrt(t_effect$mu.est.var))
effect=rbind (t_effect$mu.est,
              cred_band ,
              rep(0, length(obj$fitted_func[[cs]])) )
library(ggplot2)
library(ggtranscript)
library(biomaRt)
library(dplyr)
library(tidyr)
library(cowplot)

# Connect to Ensembl via biomaRt (current Ensembl)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve exon annotation for HSP90AA1
hsp90_exons <- getBM(
  attributes = c("chromosome_name", "exon_chrom_start", "exon_chrom_end",
                 "strand", "ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"),
  filters = "external_gene_name",
  values = "HSP90AA1",
  mart = mart
)

# Clean and format exon data
hsp90_exons <- hsp90_exons %>%
  rename(
    start = exon_chrom_start,
    end = exon_chrom_end,
    transcript = ensembl_transcript_id,
    gene = external_gene_name,
    chr = chromosome_name
  ) %>%
  mutate(
    strand = ifelse(strand == 1, "+", "-")
  )

# Step 1: Compute transcript lengths (sum of exon widths)
transcript_lengths <- hsp90_exons %>%
  group_by(transcript) %>%
  summarise(tx_length = sum(end - start + 1), .groups = "drop")

# Step 2: Join back and reorder transcript factor levels
hsp90_clean_ordered <- hsp90_exons %>%
  left_join(transcript_lengths, by = "transcript") %>%
  arrange(desc(tx_length)) %>%
  mutate(transcript = factor(transcript, levels = rev(unique(transcript)))) %>%
  arrange(transcript, start, end)

# Step 3: Generate introns per transcript
introns <- hsp90_clean_ordered %>%
  group_by(transcript) %>%
  filter(n() > 1) %>%
  to_intron()
 introns$strand="+"
# Step 4: Plot using ggtranscript with ordered transcripts
p_gene <- ggplot(hsp90_clean_ordered, aes(xstart = start, xend = end, y = transcript)) +
  geom_range(fill = "green4", height = 0.4, alpha = 0.8) +
  geom_intron(
    data = introns,strand = "-",
    aes(xstart = start, xend = end, y = transcript  ),
    arrow.min.intron.length = 250
  ) +
  labs(
    title = "",
    x = "Genomic Position",
    y = " "
  ) +
  scale_x_continuous(limits = c(102080692, 102087093), labels = scales::comma) +
  theme_cowplot() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold", angle = 0, vjust = 0.5),
    panel.grid = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

# Show the plot
print(p_gene)


positions <- seq(pos0,pos1,length.out = 1024)

markers <- obj$cs[[cs]]
j       <- which.max(obj$pip[markers])
marker  <- markers[j]
x       <- X[,marker]





read_counts <- rbind(colMeans(log1p(out$Y[x == 0,])),
                     colMeans(log1p(out$Y[x == 1,])) )
t_effect=  smashr::smash(- read_counts[1,]+read_counts[2,] , post.var=TRUE,
                         sigma= sqrt(0.05) )

library(ggplot2)

# Create data frame from the smoothed effect + band
effect_df <- data.frame(
  pos = positions,
  effect = t_effect$mu.est,
  upper = t_effect$mu.est + 1.98 * sqrt(t_effect$mu.est.var),
  lower = t_effect$mu.est - 1.98 * sqrt(t_effect$mu.est.var)
)

# Effect track plot
#p_effect <- ggplot(effect_df, aes(x = pos)) +
#  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "dodgerblue", alpha = 0.5) +
#  geom_line(aes(y = effect), color = "navyblue", linewidth = 0.7) +
#  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
#  labs(y = "Effect estimate", x = NULL) +
#  theme_cowplot() +
#  theme(
#    axis.title.x = element_blank(),
#    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
#  )

p_effect <- ggplot(effect_df, aes(x = pos)) +
  geom_line(aes(y = effect), color = "green4", linewidth = 0.7) +
  geom_line(aes(y = upper), color = "green4", linetype = "dashed", linewidth = 0.5) +
  geom_line(aes(y = lower), color ="green4", linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray30") +
  labs(y = "Effect estimate", x = NULL) +
  theme_cowplot() +
  theme(
    axis.title.x = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )


# Prepare data for ggplot
read_df <- data.frame(
  pos = rep(positions, 2),
  count = c(read_counts[1, ], read_counts[2, ]),
  genotype = rep(c("0/0", "1/1"), each = length(positions))
)

# Data track plot
p_data <- ggplot(read_df, aes(x = pos, y = count, color = genotype)) +
  geom_point(size = 0.3) + 
  scale_color_manual(values = c("darkgreen","lightgreen" )) +
  labs(y = "Avg. log1p count", x = "Genomic position (bp)", color = "Genotype") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray30") +
  theme_cowplot() +
  theme(
    legend.position = "bottom",
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )
p_data


library(cowplot)

# Remove x-axis from top two plots
p_effect_clean <- p_effect +
  theme(axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

p_data_clean <- p_data +
  theme(axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())
print(p1)  # Now should work

# Build composite plot without axis.x on first 3 panels
final_core <- plot_grid(
  p1,
  p_effect_clean,
  p_data_clean,
  p_gene,
  ncol = 1,
  align = "v",
  rel_heights = c(1.2, 1.5, 1.5, 1.5)
)


# Display it 

final_with_lines= ggdraw(final_core) +
  geom_segment(aes(x = .49, xend = 0.09,  
                   
                   y = 0.85, yend = 0.78),
               linetype = "dashed",
               color = "gray30"  )+
  geom_segment(aes(x = .575, xend = 0.96, 
                   linetype = "dashed",
                   
                   y = 0.85, yend = 0.78),
               linetype = "dashed",
               color = "gray30"  )
final_with_lines

folder_path=  paste0(getwd(),
                     "/plot/HSP90AA1/"
)
file_path <- file.path(folder_path, "HSP90AA1_cs1.pdf")
pdf(file_path, width = 8.27, height =8.27 )  # A4 in inches

final_with_lines
dev.off()
