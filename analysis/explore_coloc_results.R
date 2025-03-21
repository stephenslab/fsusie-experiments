# Some code I wrote to identify molecular trait SNPs (methylation or
# H3K27ac SNPs) that overlap with Alzheimer's disease SNPs. In the
# end, we identified two promising candidates near genes CD2AP and
# CASS4.
#
# Note: some of the files required for this analysis can be downloaded
# from https://uchicago.box.com/s/tt1vgg7vqayfthg0vsbl8gw0sdi0f1uo
#
library(data.table)

# This function is used to load the fSuSiE fine-mapping results.
read_fsusie_results <- function (filename, n = 7) {
  out <- fread(filename,sep = "\t",stringsAsFactors = FALSE,header = TRUE)
  class(out) <- "data.frame"
  out <- transform(out,chr = factor(chr))
  cols <- seq(n + 1,ncol(out))
  for (i in cols)
    out[[i]] <- factor(out[[i]])
  return(out)
}

# Read the fSuSiE results on methylation.
methyl_snps_file <- "../outputs/ROSMAP_mQTL_cs_snp_annotation.tsv.gz"
methyl_snps <- read_fsusie_results(methyl_snps_file)
methyl_snps <- methyl_snps[1:7]
methyl_snps <- transform(methyl_snps,
                         cs     = factor(cs),
                         region = factor(region),
                         study  = factor(study))

# Read the fSuSiE results on H3K27ac.
ha_snps_file <- "../outputs/ROSMAP_haQTL_cs_snp_annotation.tsv.gz"
ha_snps <- read_fsusie_results(ha_snps_file)
ha_snps <- ha_snps[1:7]
ha_snps <- transform(ha_snps,
                     cs     = factor(cs),
                     region = factor(region),
                     study  = factor(study))

# Read the coloc results.
dat <- fread("../outputs/fsusie_AD.pair_coloc.both_full_and_intersect_snps.tsv.gz",
             sep = "\t",header = TRUE,stringsAsFactors = FALSE)
class(dat) <- "data.frame"
non_ad_studies <- c("ROSMAP_DLPFC_haQTL","ROSMAP_DLPFC_mQTL")
dat <- subset(dat,
              is.element(non_ad,non_ad_studies) &
              type == "intersect")
dat_methyl <- subset(dat,non_ad == "ROSMAP_DLPFC_mQTL")
dat_ha     <- subset(dat,non_ad == "ROSMAP_DLPFC_haQTL")

# Narrow down to the cases where the evidence for colocalization is
# quite strong.
dat_methyl <- subset(dat_methyl,PP.H4.abf > 0.5)
dat_ha     <- subset(dat_ha,PP.H4.abf > 0.5)

# Narrow down to cases where the methylation CS has 4 SNPs or fewer,
# and "hit1" or "hit2" overlap with one of the SNPs in the CS.
cs_size     <- table(methyl_snps$cs)
small_cs    <- names(which(cs_size < 5))
methyl_snps <- subset(methyl_snps,is.element(cs,small_cs))
dat_methyl  <- subset(dat_methyl,
                      is.element(hit1,methyl_snps$variant_id) |
                      is.element(hit2,methyl_snps$variant_id))

# Narrow down to cases where the H3K27ac CS has 4 SNPs or fewer,
# and "hit1" or "hit2" overlap with one of the SNPs in the CS.
cs_size  <- table(ha_snps$cs)
small_cs <- names(which(cs_size < 5))
ha_snps  <- subset(ha_snps,is.element(cs,small_cs))
dat_ha   <- subset(dat_ha,
                   is.element(hit1,ha_snps$variant_id) |
                   is.element(hit2,ha_snps$variant_id))

# Display the results that remain.
options(width = 120)
cols <- c("ad","non_ad_region","hit2","PP.H4.abf","idx2")
print(dat_ha[cols])
print(dat_methyl[cols])

#
#                      ad          non_ad_region              hit2 PP.H4.abf idx2
# AD_Bellenguez_EADB_2022 chr5_82637805_88412930 chr5:87121196:C:T    0.6652    1
# AD_Bellenguez_EADB_2022 chr5_85967320_89904257 chr5:87121196:C:T    0.6653    1
#
#                    ad             non_ad_region                hit2 PP.H4.abf idx2
# AD_Kunkle_Stage1_2019    chr6_44880827_48309905   chr6:47472829:C:A    0.8582   13 *
#  AD_Wightman_..._2021    chr6_44880827_48309905   chr6:47472829:C:A    0.5194   13 *
#    AD_Bellenguez_2022    chr6_44880827_48309905   chr6:47472829:C:A    0.8957   13 *
# AD_Wightman_Full_2021 chr12_111405189_114438276 chr12:113291695:G:C    0.7752    3
# AD_Wightman_Full_2021   chr16_24031743_30613717  chr16:30001745:A:G    0.7311   13
#    AD_Bellenguez_2022   chr17_58565557_62392838  chr17:59152355:G:A    0.7412   11
#        AD_Jansen_2021   chr20_53859688_57519449  chr20:56409008:G:C    0.9058   18
#    AD_Bellenguez_2022   chr21_24561908_27573286  chr21:25870258:A:T    0.8984    7
#

# Inspect the fSuSiE results on methylation for some of these
# candidates that emerged from this analysis.
interesting_cs <- c("ROSMAP_DLPFC_mQTL:chr6_44880827_48309905:13",
                    "ROSMAP_DLPFC_mQTL:chr12_111405189_114438276:3",
                    "ROSMAP_DLPFC_mQTL:chr16_24031743_30613717:13",
                    "ROSMAP_DLPFC_mQTL:chr17_58565557_62392838:11",
                    "ROSMAP_DLPFC_mQTL:chr20_53859688_57519449:18",
                    "ROSMAP_DLPFC_mQTL:chr21_24561908_27573286:7")
cols <- c("chr","pos","cs","pip","variant_id")
print(subset(methyl_snps,
             is.element(cs,interesting_cs) &
             pip > 0.5)[cols])
#
#   chr       pos                          cs    pip          variant_id rsID       gene
#  chr6  47472829   chr6_44880827_48309905:13 0.9645   chr6:47472829:C:A rs9369695  CD2AP *
# chr12 113291695 chr12_111405189_114438276:3 1.0000 chr12:113291695:G:C rs73192839 ?
# chr16  30001745  chr16_24031743_30613717:13 0.9484  chr16:30001745:A:G rs12931955 ?
# chr17  59152355  chr17_58565557_62392838:11 1.0000  chr17:59152355:G:A rs7502947  ? 
# chr20  56409008  chr20_53859688_57519449:18 0.9999  chr20:56409008:G:C rs1884913  CASS4 *
# chr21  25870258   chr21_24561908_27573286:7 0.9247  chr21:25870258:A:T rs56117618 APP
#
# NOTE: rs1884913 is on chr 20 at 54,984,064 bp in GRCh37/hg19.
#
