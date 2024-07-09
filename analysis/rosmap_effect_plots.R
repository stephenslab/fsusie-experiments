library(fsusieR)
library(data.table)
library(ggplot2)
library(cowplot)
source("../code/rosmap_functions.R")

# Load the list of candidate SNPs.
highconf_overlap <- read.csv("../outputs/rosmap_highconf_overlap.csv",
                             stringsAsFactors = FALSE)
highconf_overlap <- transform(highconf_overlap,chr = factor(chr))

# Load the gene annotations.
gene_file <-
  file.path("../data/genome_annotations",
    "Homo_sapiens.GRCh38.103.chr.reformatted.collapse_only.gene.gtf.gz")
genes <- get_gene_annotations(gene_file)

# Load the ENCODE cCRE ("") annotations.
encode_ccre_file <- "../data/genome_annotations/encodeCcreCombined.bed.gz"
encode_ccre <- get_encode_ccre_annotations(encode_ccre_file)

# Candidate regions identified:
# 
#                 id   chr       pos     region_eqtl region_haqtl region_mqtl
# chr1:170663243:A:G  chr1 170663243 ENSG00000116132      TADB_75     TADB_75
# chr17:63440644:T:C chr17  63440644 ENSG00000008283    TADB_1205   TADB_1205
# chr19:46435226:C:T chr19  46435226 ENSG00000169515    TADB_1261   TADB_1261
#  chr3:44721338:A:G  chr3  44721338 ENSG00000163808     TADB_258    TADB_258
#
# NOTE: Start with chr3:44721338:A:G.
#

# Load the susie result.
susie <- readRDS(paste0("../data/analysis_result/Fungen_xQTL.ENSG00000163808.",
                      "cis_results_db.export.rds"))
susie <- susie[[1]][["DLPFC_DeJager_eQTL"]]

# Load the fsusie haQTL result.
fsusie_haqtl <-
   readRDS(paste0("../data/analysis_result/ROSMAP_haQTL.chr3_43915257_",
                  "48413435.fsusie_mixture_normal_top_pc_weights.rds"))
fsusie_haqtl <- fsusie_haqtl[[1]][[1]]

# Load the fsusie mQTL result.
fsusie_mqtl <-
  readRDS(paste0("../data/analysis_result/ROSMAP_mQTL.chr3_43915257_",
                 "48413435.fsusie_mixture_normal_top_pc_weights.rds"))
fsusie_mqtl <- fsusie_mqtl[[1]][[1]]

# Create the annotated effect plot.
chr <- "chr3"
mb0 <- 44.25
mb1 <- 44.9
snp_pos <- 44721338/1e6
fsusie_haqtl$fsusie_result$outing_grid <-
  fsusie_haqtl$fsusie_result$outing_grid/1e6
fsusie_mqtl$fsusie_result$outing_grid <-
  fsusie_mqtl$fsusie_result$outing_grid/1e6
p1 <- plot_susiF_effect(fsusie_haqtl$fsusie_result,effect = 1,
                        title = "H3K27ac") +
  geom_vline(xintercept = snp_pos,color = "gray",linetype = "dotted") +
  xlim(mb0,mb1) +
  xlab("") 
p2 <- plot_susiF_effect(fsusie_mqtl$fsusie_result,effect = 8,
                        title = "methylation") +
  geom_vline(xintercept = snp_pos,color = "gray",linetype = "dotted") +
  xlim(mb0,mb1) +
  xlab("")
gene_dat <- subset(genes,chromosome == chr &
                   start < mb1*1e6 &
                   end > mb0*1e6)
gene_dat <- transform(gene_dat,
                      start = start/1e6,
                      end = end/1e6)
rows <- which(gene_dat$strand == "-")
gene_start <- gene_dat[rows,"start"]
gene_end   <- gene_dat[rows,"end"]
gene_dat[rows,"start"] <- gene_end
gene_dat[rows,"end"]   <- gene_start
p3 <- ggplot(gene_dat,aes(x = start,xend = end,y = gene_name)) +
  geom_segment(linewidth = 0.35,arrow = arrow(length = unit(0.1,"cm"))) +
  labs(x = "",y = "") +
  geom_vline(xintercept = snp_pos,color = "gray",linetype = "dotted") +
  theme_cowplot(font_size = 8) 
ucsc_colors <- c("dodgerblue","gold","darkorange","pink","red")
annotation_dat <- subset(encode_ccre,
                         chrom == chr &
                         chromStart < mb1*1e6 &
                         chromEnd > mb0*1e6)
annotation_dat <- transform(annotation_dat,
                            chromPos = (chromStart + chromEnd)/2)
annotation_dat <- transform(annotation_dat,
                            chromStart = chromStart/1e6,
                            chromEnd   = chromEnd/1e6,
                            chromPos   = chromPos/1e6)
p4 <- ggplot(annotation_dat,aes(x = chromPos,y = ucscLabel,
                                color = ucscLabel)) +
  geom_point(shape = 17,size = 0.85) +
  geom_vline(xintercept = snp_pos,color = "gray",linetype = "dotted") +
  scale_color_manual(values = ucsc_colors) +
  labs(x = "base-pair position (Mb)") +
  theme_cowplot(font_size = 9)
plot_grid(p1,p2,p3,p4,nrow = 4,ncol = 1,align = "v",axis = "lr")
