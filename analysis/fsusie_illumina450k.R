# Short script to take a quick look at the fSuSiE mSNPs that overlap
# with the CpG probes.
#
# Downloaded IlluminaHumanMethylation450K.rda from
# https://github.com/Yang9704/MethylCallR
#
# Downloaded 12864_2013_7006_MOESM2_ESM.csv.gz from
# https://doi.org/10.1186/1471-2164-15-51
#
library(data.table)
library(ggplot2)
library(cowplot)

# Import the Illumina HumanMethylation450 CpG probe information.
load("../data/IlluminaHumanMethylation450K.rda")
probes <- IlluminaHumanMethylation450K
probes <- subset(probes,!is.na(Start_hg38) & !is.na(End_hg38))
probes <- subset(probes,End_hg38 > Start_hg38)
probes <- probes[c("CHR_hg38","Start_hg38","End_hg38")]
names(probes) <- c("chr","start","end")
probes <- transform(probes,width = end - start)
probes$cpg <- rownames(probes)
rownames(probes) <- NULL

# Import the Naeem et al recommendations for CpGs to keep and discard.
naeem <- read.csv("../data/12864_2013_7006_MOESM2_ESM.csv.gz",
                  stringsAsFactors = FALSE)

# Load the fSuSiE mSNP results.
read_enrichment_results <- function (filename, n) {
  out <- fread(filename,sep = "\t",stringsAsFactors = FALSE,header = TRUE)
  class(out) <- "data.frame"
  out <- transform(out,chr = factor(chr))
  if (ncol(out) > n) {
    cols <- seq(n + 1,ncol(out))
    for (i in cols)
      out[[i]] <- factor(out[[i]])
  }
  return(out)
}
methyl_snps_fsusie_file <- "../outputs/ROSMAP_mQTL_cs_snp_annotation.tsv.gz"
methyl_snps_fsusie <- read_enrichment_results(methyl_snps_fsusie_file,n = 7)
methyl_snps_fsusie <- transform(methyl_snps_fsusie,
                                cs     = factor(cs),
				region = factor(region),
				study  = factor(study))
methyl_snps_fsusie <- methyl_snps_fsusie[1:7]
methyl_snps_fsusie <- transform(methyl_snps_fsusie,
                                id = paste(chr,pos,sep = "_"))

# Flag the SNPs that are inside the CpG probes.
get_probe_sites <- function (probes)
  c(with(subset(probes,width == 1),
         c(paste(chr,start,sep = "_"),
           paste(chr,start + 1,sep = "_"))),
    with(subset(probes,width == 2),
         c(paste(chr,start,sep = "_"),
           paste(chr,start + 1,sep = "_"),
           paste(chr,start + 2,sep = "_"))))
cpgs_keep    <- subset(naeem,Flag.discard.keep. == "keep")$probe
cpgs_discard <- subset(naeem,Flag.discard.keep. == "discard")$probe
probe_sites_keep <- get_probe_sites(subset(probes,is.element(cpg,cpgs_keep)))
probe_sites_discard <- get_probe_sites(subset(probes,
                                              is.element(cpg,cpgs_discard)))
methyl_snps_fsusie <-
  transform(methyl_snps_fsusie,
            in_probe_keep    = is.element(id,probe_sites_keep),
            in_probe_discard = is.element(id,probe_sites_discard))
length(unique(subset(methyl_snps_fsusie,in_probe_keep)$id))
# 57
length(unique(subset(methyl_snps_fsusie,in_probe_discard)$id))
# 942

# Flag the SNPs from the association tests that are inside the CpG probes.
methyl_cpg_assoc_file <- "../outputs/ROSMAP_mQTL_qtl_snp_qval0.05.tsv.gz"
methyl_cpg_assoc <- read_enrichment_results(methyl_cpg_assoc_file,n = 8)
ids <- with(methyl_cpg_assoc,sprintf("%s_%d",chr,pos))
ids <- unique(ids)
summary(is.element(ids,probe_sites_keep))
#     FALSE     TRUE
#  10533116     1463
summary(is.element(ids,probe_sites_discard))
#    FALSE     TRUE
# 10520840    13739
 
# "Blacklist" the fSuSiE CSs.
methyl_snps_fsusie <- transform(methyl_snps_fsusie,
                                in_probe = in_probe_keep | in_probe_discard)
blacklisted <- with(methyl_snps_fsusie,
                    tapply(in_probe,cs,any))
blacklisted <- names(which(blacklisted))
methyl_snps_fsusie <- transform(methyl_snps_fsusie,
                                blacklisted = is.element(cs,blacklisted))

# Look at the distribution of the CS sizes for CSs overlapping and not
# overlapping a CpG probe.
bins <- c(0,1,2,5,10,20,Inf)
cs_sizes_out <- table(factor(subset(methyl_snps_fsusie,!blacklisted)$cs))
cs_sizes_in  <- table(factor(subset(methyl_snps_fsusie,blacklisted)$cs))
cs_sizes_out <- cut(as.numeric(cs_sizes_out),bins)
cs_sizes_in  <- cut(as.numeric(cs_sizes_in),bins)
levels(cs_sizes_in) <- bins[-1]
levels(cs_sizes_out) <- bins[-1]
p1 <- ggplot(data.frame(cs_size = cs_sizes_out),aes(x = cs_size)) +
  geom_histogram(stat = "count",color = "white",fill = "darkblue",
                 width = 0.65) +
  labs(x = "CS size",y = "number of CSs",
       title = "not overlapping a CpG probe") +
  theme_cowplot(font_size = 10) +
  theme(plot.title = element_text(size = 10,face = "plain"))
p2 <- ggplot(data.frame(cs_size = cs_sizes_in),aes(x = cs_size)) +
  geom_histogram(stat = "count",color = "white",fill = "darkblue",
                 width = 0.65) +
  labs(x = "CS size",y = "number of CSs",title = "overlapping a CpG probe") +
  theme_cowplot(font_size = 10) +
  theme(plot.title = element_text(size = 10,face = "plain"))
plot_grid(p1,p2,nrow = 1,ncol = 2)
