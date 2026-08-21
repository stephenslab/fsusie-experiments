# https://kevinlkx.github.io/fiber-seq-analysis/
# example_data_topic_wavelet_analysis.html
#
# Data downloaded from:
# /project/spott/kevinluo/Fiber_seq/results/example_regions/processed_data/
#
# TO DO: Convert this analysis to R Markdown.
# 
library(ggplot2)
library(cowplot)
library(susieR)
library(fsusieR)
set.seed(1)
dat <- readRDS("../data/processed_region_data_chr2_212990848_212991848.rds")
dat <- readRDS("../data/processed_region_data_chr2_212990848_212991848.rds")
# This is a function we will use in the code chunks below to visualize
# the results of the fSuSiE fine-mapping analysis.
pip_plot <- function (fit, causal_snps = NULL, xlab = "") {
  cs_colors <- c("dodgerblue","limegreen","darkorange")
  pip <- fit$pip
  n   <- length(pip)
  dat <- data.frame(pos = 1:n,pip = pip,CS = as.character(NA))
  num_cs <- length(fit$cs)
  for (i in num_cs:1) {
    j <- fit$cs[[i]]
	dat[j,"CS"] <- paste0("CS",i)
  }
  dat_cs <- subset(dat,!is.na(CS))
  out <- ggplot(dat,aes(x = pos,y = pip)) +
    geom_point() +
	geom_point(data = dat_cs,shape = 1,size = 3,
	           mapping = aes(color = CS)) +
	scale_x_continuous(breaks = NULL) +
	scale_color_manual(values = cs_colors) + 
	labs(x = xlab,y = "PIP") +
    theme_cowplot(font_size = 12)
  if (!is.null(causal_snps)) {
	dat_causal <- dat[causal_snps,]
    out <- out + geom_point(data = dat_causal,color = "red")
  }
  return(out)
}

# Get the SNP genotypes.
geno <- dat$genotype_dosage_mat
maf  <- colMeans(geno)/2
maf  <- pmin(maf,1 - maf)
j    <- which(maf > 0.05)
geno <- geno[,j]

# Get the pseudobulked read counts.
pheno <- NULL
metA <- as.matrix(dat$metA_mat)
rows <- which(!apply(is.na(metA),1,any))
metA <- metA[rows,]
sample_ids <- dat$rids_df$sample_ID[rows]
for (sample_id in levels(sample_ids)) {
  i <- which(sample_ids == sample_id)
  pheno <- rbind(pheno,colMeans(metA[i,]))
}
rownames(pheno) <- levels(sample_ids)
print(all(rownames(geno) == rownames(pheno)))

fit <- susiF(pheno,geno,L = 2,quantile_trans = FALSE,
             filter_cs = FALSE,maxit = 100,
             post_processing = "smash")

pip_plot(fit)
fit$cs
R <- cor(geno)
j <- fit$cs[[1]]
R[j,j]
j <- fit$cs[[2]]
R[j,j]

# Plot the estimated effects of the two causal SNPs.
x <- as.numeric(colnames(pheno))
out1 <- get_fitted_effect(fit,l = 1,cred_band = TRUE,alpha = 0.05)
out2 <- get_fitted_effect(fit,l = 2,cred_band = TRUE,alpha = 0.05)
plot(x,out1$effect,type = "l",col = "dodgerblue",xlab = "base-pair position",
     ylab = "effect")
lines(x,out2$effect,col = "orangered")

j <- fit$cs[[1]][1]
i0 <- which(geno[,j] == 0)
i1 <- which(geno[,j] == 1)
i2 <- which(geno[,j] == 2)
plot(x,colMeans(pheno[i0,]),type = "l",col = "darkblue",
     xlab = "base-pair position")
lines(x,colMeans(pheno[i1,]),type = "l",col = "cyan")
lines(x,colMeans(pheno[i2,]),type = "l",col = "orangered")

# TO DO NEXT: Try analyzing 
# processed_region_data_chr1_37690171_37691171.rds.
