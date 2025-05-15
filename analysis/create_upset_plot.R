# TO DO: Explain here what this script does, and how to use it.
library(ggplot2)
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

# Create the upset plot.
# TO DO.
