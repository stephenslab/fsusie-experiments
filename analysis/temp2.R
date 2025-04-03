# x <- factor(methyl_cpg_assoc$variant_id)
# cpgs_per_snp_assoc <- tapply(methyl_cpg_assoc$molecular_trait_id,x,
#                              function (x) length(unique(x)))
# rm(x)
nodup_cs <- names(which(table(methyl_snps_fsusie_nodup$cs) > 0))
methyl_cpg_fsusie_nodup <- subset(methyl_cpg_fsusie,is.element(cs,nodup_cs))
cpgs_per_snp_fsusie <-
  with(methyl_cpg_fsusie_nodup,
       tapply(ID,cs,function (x) length(unique(x))))
pdat1 <- data.frame(x = cpgs_per_snp_assoc)
pdat2 <- data.frame(x = cpgs_per_snp_fsusie)
pdat1 <- subset(pdat1,x <= 25)
pdat2 <- subset(pdat2,x <= 250)
p1 <- ggplot(pdat1,aes(x = x)) +
  geom_histogram(color = "white",fill = "darkblue",bins = 25) +
  labs(x = "number of CpGs",y = "number of SNPs",
       title = "SNP-CpG association tests") +
  theme_cowplot(font_size = 10) +
  theme(plot.title = element_text(size = 10,face = "plain"))
p2 <- ggplot(pdat2,aes(x = x)) +
  geom_histogram(color = "white",fill = "darkblue",bins = 25) +
  labs(x = "number of CpGs",y = "number of CSs",title = "fSuSiE") +
  theme_cowplot(font_size = 10) +
  theme(plot.title = element_text(size = 10,face = "plain"))
print(plot_grid(p1,p2,nrow = 1,ncol = 2))
