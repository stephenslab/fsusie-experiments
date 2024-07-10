library(biomaRt)
highconf_overlap <- read.csv("../outputs/rosmap_highconf_overlap.csv",
                             stringsAsFactors = FALSE)

# Note that this is hg38.
snpMart <- useEnsembl(biomart = "snps",dataset = "hsapiens_snp",
                      mirror = "www")
dat <- data.frame(CHR   = as.numeric(substr(highconf_overlap$chr,4,5)),
                  START = highconf_overlap$pos,
                  END   = highconf_overlap$pos)
coords <- apply(dat,1,paste,collapse = ":")

# Submit the query.
res <- getBM(attributes = c("refsnp_id","chr_name","chrom_start","chrom_end",
                            "allele"),
             filters = "chromosomal_region",
             values = coords[1:10],
             mart = snpMart)  
