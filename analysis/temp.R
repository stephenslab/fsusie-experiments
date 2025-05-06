temp <- subset(methyl_cs1snp_fsusie,maf > 0.03)
methyl_cpg_fsusie_cs1snp <- subset(methyl_cpg_fsusie,
                                   is.element(cs,temp$cs))
rows <- match(methyl_cpg_fsusie_cs1snp$cs,temp$cs)
methyl_cpg_fsusie_cs1snp$variant_pos <- temp[rows,"pos"]
methyl_cpg_fsusie_cs1snp <- transform(methyl_cpg_fsusie_cs1snp,
                                      cs = factor(cs),
                                      dist = peak_start - variant_pos)
pdat <- tapply(methyl_cpg_fsusie_cs1snp$dist,methyl_cpg_fsusie_cs1snp$cs,
               function (x) {
                 i <- which.min(abs(x))
                 return(x[i])
               })
pdat <- data.frame(dist = pdat)
pdat <- transform(pdat,
  dist = cut(dist,c(-Inf,-1e6,-1e5,-1e4,-1e3,-100,
                    100,1e3,1e4,1e5,1e6,+Inf)))
p <- ggplot(pdat,aes(x = dist)) +
  geom_bar(color = "white",fill = "darkblue",bins = 64,width = 0.5) +
  labs(x = "CpG - mSNP distance (base-pairs)",
       y = "number of mSNPs") +
  theme_cowplot(font_size = 10) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) 
print(p)
