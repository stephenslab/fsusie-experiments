# Get the MAFs for all SNPs with MAF > 0.1%, and save as an .RData
# file.
library(tools)
library(data.table)
afreq_file <- "../data/ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.afreq.gz"
system(paste("gzcat",afreq_file,"| cut -f 1,6 > temp1.txt"))
system(paste("gzcat",afreq_file,"| cut -f 2 | cut -d ':' -f 2 |",
             "cut -d '_' -f 1 > temp2.txt"))
system(paste("paste temp1.txt temp2.txt > rosmap.afreq"))			 
afreq <- fread("rosmap.afreq",sep = "\t",stringsAsFactors = FALSE,
               header = TRUE,nrows = Inf)
class(afreq) <- "data.frame"
names(afreq) <- c("chr","maf","pos")
afreq <- transform(afreq,
                   chr = factor(chr),
				   maf = pmin(maf,1 - maf))
afreq <- subset(afreq,maf > 0.001)
levels(afreq$chr) <- paste0("chr",1:22)
afreq$id <- sprintf("%s_%d",as.character(afreq$chr),afreq$pos)
save(list = "afreq",file = "afreq.RData")
resaveRdaFiles("afreq.RData")
