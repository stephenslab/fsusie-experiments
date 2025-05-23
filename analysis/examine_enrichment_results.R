library(ggplot2)
library(cowplot)
# trait <- "ROSMAP_mQTL"
trait <- "ROSMAP_haQTL"
infile <- paste("../outputs/all.fsusie.susie_toppc1.haQTL.mQTL.enrichment",
                "results_encode.tsv.gz",sep = "_")
dat <- read.table(infile,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
dat <- dat[c("Annotation","Enrichment_log2","Enrichment_SE_log2",
             "Enrichment_P_value","context","method","resource")]
dat <- transform(dat,
                 context  = factor(context), 
                 method   = factor(method),
                 resource = factor(resource))
dat <- subset(dat,context == trait)             
annotations <- c("CTCF_Hoffman",
                 "Coding_UCSC",
                 # "Conserved_LindbladToh",
                 "DHS_peaks_Trynka",
                 # "DNaseI",
                 "Enhancer_Andersson",
                 "Enhancer_Hoffman",
                 # "FetalDHS_Trynka",
                 "PromoterFlanking_Hoffman",
                 "Promoter_UCSC",
                 "Repressed_Hoffman",
                 "SuperEnhancer_Hnisz",
                 "TSS_Hoffman",
                 "UTR_3_UCSC",
                 "UTR_5_UCSC",
                 "WeakEnhancer_Hoffman",
                 "Intron_UCSC")
dat <- subset(dat,is.element(Annotation,annotations))
dat <- transform(dat,Annotation = factor(Annotation,annotations))
x   <- tapply(dat$Enrichment_log2,dat$Annotation,max)
annotations <- names(sort(x))
dat <- transform(dat,Annotation = factor(Annotation,annotations))
p <- ggplot(dat,aes(x = Enrichment_log2,y = Annotation,color = method,
                    xmin = Enrichment_log2 - Enrichment_SE_log2,
                    xmax = Enrichment_log2 + Enrichment_SE_log2)) +
  geom_point(shape = 20,size = 3) +
  geom_errorbarh(height = 0) +
  geom_vline(xintercept = 0,linetype = "dotted") +
  scale_color_manual(values = c("magenta","dodgerblue")) +
  labs(x = "enrichment",y = "",title = trait) +
  theme_cowplot(font_size = 10)
print(p)
outfile <- paste0("enrichment_ldsc_",trait,".pdf")
ggsave(outfile,p,height = 2.5,width = 4.5)
