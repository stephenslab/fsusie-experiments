library(data.table)
dat <- fread("../data/ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.afreq.gz",
             sep = "\t",stringsAsFactors = FALSE,header = TRUE,nrows = Inf)
class(dat) <- "data.frame"
hist(log10(dat$ALT_FREQS),n = 64)
