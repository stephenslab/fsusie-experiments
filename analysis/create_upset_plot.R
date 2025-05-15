# Create the "upset plot" showing the overlap among the different
# molecular trait SNPs (RNA expression, protein expression, H3K27ac,
# methylation).
library(ggplot2)
library(tibble)
library(cowplot)
library(ggupset)

# Get the fSuSiE 1-SNP CSs for H3K27ac and methylation.
ha_fsusie     <- readRDS("../outputs/h3k27ac_cs1snp_fsusie.rds")
methyl_fsusie <- readRDS("../outputs/methyl_cs1snp_fsusie.rds")

# Get the SuSiE 1-SNP CSs for RNA expression.
rnaseq_susie  <- readRDS("../outputs/rnaseq_susie.rds")
gene_cs <- with(rnaseq_susie,paste(gene,cs_coverage_0.95_min_corr,sep = "_"))
gene_cs <- factor(gene_cs)
cs_size <- table(gene_cs)
cs1snp  <- names(which(cs_size == 1))
rows    <- which(is.element(gene_cs,cs1snp))
rnaseq_susie <- rnaseq_susie[rows,]

# Get the SuSiE 1-SNP CSs for protein expression.
protein_susie <- readRDS("../outputs/protein_susie.rds")
gene_cs <- with(protein_susie,
                paste(gene,cs_coverage_0.95_min_corr,sep = "_"))
gene_cs <- factor(gene_cs)
cs_size <- table(gene_cs)
cs1snp  <- names(which(cs_size == 1))
rows    <- which(is.element(gene_cs,cs1snp))
protein_susie <- protein_susie[rows,]

# Add a "chr_pos" column to each data frame.
rnaseq_susie  <- transform(rnaseq_susie,chr_pos = paste(chr,pos,sep = "_"))
protein_susie <- transform(protein_susie,chr_pos = paste(chr,pos,sep = "_"))
ha_fsusie     <- transform(ha_fsusie,chr_pos = paste(chr,pos,sep = "_"))
methyl_fsusie <- transform(methyl_fsusie,chr_pos = paste(chr,pos,sep = "_"))

# Create a "cs" column for the SuSiE results.
rnaseq_susie <- 
  transform(rnaseq_susie,cs = paste(gene,cs_coverage_0.95_min_corr,sep="_"))
protein_susie <- 
  transform(protein_susie,cs = paste(gene,cs_coverage_0.95_min_corr,sep="_"))

# Keep only the columns that are needed to create the upset plot.
rnaseq_susie  <- rnaseq_susie[c("gene","chr_pos")]
protein_susie <- protein_susie[c("gene","chr_pos")]
ha_fsusie     <- ha_fsusie[c("cs","chr_pos")]
methyl_fsusie <- methyl_fsusie[c("cs","chr_pos")]
names(rnaseq_susie)  <- c("rnaseq_cs","chr_pos")
names(protein_susie) <- c("protein_cs","chr_pos")
names(ha_fsusie)     <- c("ha_cs","chr_pos")
names(methyl_fsusie) <- c("methyl_cs","chr_pos")

# Combine all CSs into a single data frame.
pdat <- merge(rnaseq_susie,protein_susie,all = TRUE,sort = FALSE)
pdat <- merge(pdat,ha_fsusie,all = TRUE,sort = FALSE)
pdat <- merge(pdat,methyl_fsusie,all = TRUE,sort = FALSE)

# Create a "list" column containing the affected molecular traits.
n <- nrow(pdat)
traits <- vector("list",n)
trait_names <- c("rnaseq_cs","protein_cs","ha_cs","methyl_cs")
for (i in 1:n) {
  x <- !is.na(pdat[i,trait_names])
  traits[[i]] <- trait_names[which(x)]
}
pdat <- as_tibble(pdat)
pdat$traits <- traits

# Create the upset plot.
p <- ggplot(pdat,aes(x = traits)) +
  geom_bar(color = "black",fill = "black",width = 0.5) +
  geom_text(stat = "count",aes(label = after_stat(count)),vjust = -1,
            size = 3) +
  scale_x_upset(order_by = "degree",sets = trait_names) +
  scale_y_continuous(breaks = seq(0,10000,1000),limits = c(0,7000)) +
  labs(x = "",y = "number of 1-SNP CSs") +
  theme_cowplot(font_size = 10)
ggsave("upset_plot.pdf",p,height = 4,width = 4)
