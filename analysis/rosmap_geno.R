# First run
# gzcat ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.afreq.gz \
#   | cut -f 1,6 > temp.tsv
library(data.table)
library(ggplot2)
library(cowplot)
dat <- fread("../data/temp.tsv",sep = "\t",stringsAsFactors = FALSE,
             header = TRUE,nrows = Inf)
class(dat) <- "data.frame"
x <- dat$ALT_FREQS
x <- cut(x,c(0,0.0001,0.001,0.005,0.01,0.05,Inf))
counts <- table(x)
pdat <- data.frame(count = as.vector(counts),
                   bins  = names(counts))
print(ggplot(pdat,aes(x = bins,y = count)) +
      geom_col(width = 0.4,color = "white",fill = "royalblue") +
      labs(x = "MAF",y = "number of SNPs") +
      theme_cowplot(font_size = 9) +
      theme(axis.text.x = element_text(angle = 45,hjust = 1)))
