library(ggplot2)
library(cowplot)
trait <- "mQTL"
# trait <- "haqtl"
# infile1 <-
#   file.path("../outputs",
#             paste(trait,"qtl_snp_qvalue0.05.enrichment_results_summary.tsv.gz",sep = "_"))
infile1 <-
  file.path("../outputs",
    paste(trait,"cs_snp_toppc1_pip0.5.enrichment_results_summary.tsv.gz",
          sep = "_"))
infile2 <-
  file.path("../outputs",
    paste(trait,"cs_snp_pip0.5.enrichment_results_summary.tsv.gz",sep = "_"))
dat1 <- read.table(infile1,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
dat2 <- read.table(infile2,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
hist(dat1$Enrichment,n = 32)
hist(dat2$Enrichment,n = 32)
plot(dat1$Enrichment,dat2$Enrichment,pch = 20)
abline(a = 0,b = 1,col = "magenta",lty = "dotted")
annotations <- c("Enhancer_Andersson.to_hg38",
                 "Enhancer_Hoffman.to_hg38",
                 "SuperEnhancer_Hnisz.to_hg38",
                 "WeakEnhancer_Hoffman.to_hg38",
                 "PromoterFlanking_Hoffman.to_hg38",
                 "Promoter_UCSC.to_hg38",
                 "TSS_Hoffman.to_hg38",
                 "UTR_3_UCSC.to_hg38",
                 "UTR_5_UCSC.to_hg38",
                 "Intron_UCSC.to_hg38",
                 "UTR_3_UCSC.to_hg38",
                 "UTR_5_UCSC.to_hg38",
                 "DNaseI.to_hg38",
                 "E081-DNase.macs2.narrowPeak.to_hg38",
                 "DHS_peaks_Trynka.to_hg38",
                 "FetalDHS_Trynka.to_hg38",
                 "CTCF_Hoffman.to_hg38",
                 "Repressed_Hoffman.to_hg38")
rows <- which(is.element(dat1$Annotation,annotations))
plot(dat1[rows,"Enrichment"],dat2[rows,"Enrichment"],pch = 20)
abline(a = 0,b = 1,col = "magenta",lty = "dotted")
dat1 <- dat1[rows,]
dat2 <- dat2[rows,]
pdat <- rbind(data.frame(method = "SuSiE-topPC",
                         annotation = dat1$Annotation,
                         enrichment = dat1$Enrichment,
                         stringsAsFactors = FALSE),
              data.frame(method = "fSuSiE",
                         annotation = dat2$Annotation,
                         enrichment = dat2$Enrichment,
                         stringsAsFactors = FALSE))
pdat <- transform(pdat,
                  method     = factor(method),
                  annotation = factor(annotation),
                  enrichment = log10(enrichment))
print(ggplot(pdat,aes(x = enrichment,y = annotation,color = method)) +
      geom_point(shape = 20,size = 3) +
      scale_color_manual(values = c("magenta","dodgerblue")) +
      xlim(c(-0.5,1.1)) + 
      labs(y = "",title = trait) +
      theme_cowplot(font_size = 10) +
      theme(panel.grid.major.y = element_line(linewidth = 0.5,
                                              color = "darkgray")))

