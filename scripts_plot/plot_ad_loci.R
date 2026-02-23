library(data.table)
library(ggplot2)
library(ggrepel)
library(cowplot)

# This is a script to plot the results from the 14 methylation and
# histone acetylation loci that overlap with AD loci according to
# coloc.
ad_loci <- data.frame(trait = c(rep("haQTL",4),rep("mQTL",14)),
                      region = c("chr1_205117782_208795513",
                                 "chr5_82637805_88412930",
                                 "chr5_85967320_89904257",
                                 "chr8_98960936_104073051",
                                 "chr6_44880827_48309905",
                                 "chr10_54746550_59267894",
                                 "chr12_109831778_113509918",
                                 "chr12_111405189_114438276",
                                 "chr12_111570379_117277693",
                                 "chr16_21437007_24596316",
                                 "chr16_24031743_30613717",
                                 "chr16_27341433_34991551",
                                 "chr17_58565557_62392838",
                                 "chr18_59611772_63187118",
                                 "chr18_62476957_68137751",
                                 "chr20_53859688_57519449",
                                 "chr21_24561908_27573286",
                                 "chr21_25835755_30175889"))

# Load the relevant coloc results.
coloc_res <- fread("../outputs/ROSMAP_fsusie_AD_coloc.tsv.gz",
                   sep = "\t",header = TRUE)
class(coloc_res) <- "data.frame"

# Load the relevant AD GWAS results.
ad_gwas <- fread("../outputs/fsusie_AD_overlap_sumstat.tsv.gz",
                 sep = "\t",header = TRUE)
class(ad_gwas) <- "data.frame"

# Add some columns to ad_loci for storing the coloc results.
ad_loci <- cbind(ad_loci,
                 data.frame(study       = "",
                            PP.H0.abf   = 0,
                            PP.H1.abf   = 0,
                            PP.H2.abf   = 0,
                            L_PP.H3.abf = 0,
                            L_PP.H4.abf = 0))

# Repeat for each of the selected loci.
n <- nrow(ad_loci)
for (i in 1:n) {
  trait <- ad_loci[i,"trait"]
  region <- ad_loci[i,"region"]
  cat("trait =",trait,"|","region =",region,"\n")
  
  # Extract the coloc results for the region.
  cols <- c("PP.H0.abf","PP.H1.abf","PP.H2.abf","L_PP.H3.abf","L_PP.H4.abf")
  res <- subset(coloc_res,xQTL == sprintf("ROSMAP_%s@%s",trait,region))
  ad_loci[i,"study"] <- res[1,"AD"]
  ad_loci[i,cols] <- res[1,cols]
  
  # Plot the Alzheimer's disease GWAS results.
  dat1 <- subset(ad_gwas,
                 type == ad_loci[i,"trait"] &
                 region == ad_loci[i,"region"])
  dat1 <- transform(dat1,
                    pos = pos/1e6,
                    pvalue = -log10(pvalue))
  chr <- dat1[i,"chrom"]
  study <- dat1[i,"AD"]
  p1 <- ggplot(dat1,aes(x = pos,y = pvalue,)) +
    geom_point() +
    labs(x = sprintf("base-pair position on chromosome %d (Mb)",chr),
         y = "-log10 p-value",
         title = study) +
    theme_cowplot(font_size = 10)
  print(p1)
  
  # Save the plots in a PDF.
  # TO DO.
}
