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
                                 "chr21_25835755_30175889"),
                      pos0 = c(207,rep(0,17)),
                      pos1 = c(208,rep(0,17)))

# Load the allele frequency data.
load("../data/afreq.RData")
afreq <- subset(afreq,maf >= 0.05)

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
  
  # Plot the Alzheimer's disease GWAS results and the fsusie
  # fine-mapping results (PIPs and CSs). Note that the CSs are
  # filtered so that the sentinel SNP always has an MAF > 5%.
  dat1 <- subset(ad_gwas,
                 type == ad_loci[i,"trait"] &
                 region == ad_loci[i,"region"])
  dat1$label = ""
  j <- which.min(dat1$pvalue)
  dat1[j,"label"] <- prettyNum(dat1[j,"pos"],big.mark = ",",scientific = FALSE)
  dat1 <- transform(dat1,
                    pos    = pos/1e6,
                    pvalue = -log10(pvalue))
  chromosome <- dat1[i,"chrom"]
  study <- dat1[i,"AD"]
  
  rdsname <-
    file.path("../outputs/fsusie_ad_loci",
              sprintf("ROSMAP_%s.%s.fsusie_mixture_normal_top_pc_weights.rds",
                      trait,region))
  fsusie <- readRDS(rdsname)$fsusie_obj
  pip    <- fsusie$pip
  cs     <- fsusie$sets$cs
  dat2   <- data.frame(pip   = pip,
                       pos   = as.numeric(sapply(names(pip),
                                 function (x) unlist(strsplit(x,":"))[2])),
                       cs    = 0,
                       label = "",
                       variant_id = names(pip))
  rownames(dat2) <- names(pip)
  for (j in 1:length(cs)) {
    snps <- names(cs[[j]])
    sentinel_snp <- names(which.max(pip[snps]))
    if (sum(with(afreq,chr == paste0("chr",chromosome) &
                       pos == dat2[sentinel_snp,"pos"]) > 0)) {
      dat2[snps,"cs"] <- j
      dat2[sentinel_snp,"label"] <-
        prettyNum(dat2[sentinel_snp,"pos"],big.mark = ",",scientific = FALSE)
    } else {
      snps <- setdiff(rownames(dat2),snps)
      dat2 <- dat2[snps,]
    }
  }
  dat2 <- transform(dat2,cs = factor(cs))
  dat2 <- transform(dat2,pos = pos/1e6)
  cs_labels <- sprintf("CS%d (%d SNPs)",1:nlevels(dat2$cs) - 1,table(dat2$cs))
  j <- which(table(dat2$cs) == 1)
  cs_labels[j] <- sprintf("CS%d (1 SNP)",j - 1)
  cs_labels[1] <- "none"
  levels(dat2$cs) <- cs_labels
  rownames(dat2) <- NULL
  dat1$cs <- "none"
  for (j in cs_labels) {
    snps <- subset(dat2,cs == j)$variant_id
    snps <- substr(snps,4,100)
    rows <- which(is.element(dat1$variant_id,snps))
    dat1[rows,"cs"] <- j
  }
  dat1 <- transform(dat1,cs = factor(cs,cs_labels))
  levels(dat1$cs) <- substr(levels(dat1$cs),1,4)
  dat1 <- dat1[order(dat1$cs),]
  cs_colors <- c("black","red","dodgerblue","limegreen","darkorange",
                 "gold","peru","violet","#FFB300","#803E75","#FF6800",
                 "#A6BDD7","#C10020","#CEA262","#817066","#007D34")
  p1 <- ggplot(dat1,aes(x = pos,y = pvalue,label = label,color = cs)) +
      geom_point(show.legend = TRUE) +
    geom_text_repel(size = 2.5,color = "black",segment.color = "black",
                    min.segment.length = 0,max.overlaps = Inf) +
    scale_x_continuous(breaks = seq(floor(min(dat1$pos)),
                                    ceiling(max(dat1$pos)),0.5)) +
    scale_color_manual(values = cs_colors,drop = FALSE) +
    labs(x = sprintf("base-pair position on chromosome %d (Mb)",chromosome),
         y = "-log10 p-value",color = "CS",
         title = paste0(study,", PP(H4) = ",
                        round(coloc_res[i,"L_PP.H4.abf"],digits = 3))) +
    theme_cowplot(font_size = 9) +
    theme(plot.title = element_text(size = 9,face = "plain"))
  p2 <- ggplot(dat2,aes(x = pos,y = pip,color = cs,label = label)) +
    geom_point(show.legend = TRUE) +
    geom_text_repel(size = 2.5,color = "black",segment.color = "black",
                    min.segment.length = 0,max.overlaps = Inf) +
    scale_x_continuous(breaks = seq(min(dat2$pos),max(dat2$pos),0.5)) +
    scale_y_continuous(limits = c(0,1.1),breaks = c(0,0.5,1)) +
    scale_color_manual(values = cs_colors,drop = FALSE) +
    labs(x = sprintf("base-pair position on chromosome %d (Mb)",chromosome),
         y = "PIP",color = "CS",
         title = ifelse(trait == "haQTL","haSNPs in DFPLC","mSNPs in DFPLC")) +
    theme_cowplot(font_size = 9) +
    theme(plot.title = element_text(size = 9,face = "plain"))
    
  # Save the plots to a PDF.
  # print(plot_grid(p1,p2,nrow = 2,ncol = 1,align = "v"))
  pdfname <- sprintf("%s_%s_plots.pdf",trait,region)
  ggsave(pdfname,
         plot_grid(p1,p2,nrow = 2,ncol = 1,align = "v"),
         height = 3,width = 5)
}
