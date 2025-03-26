library(ggplot2)
library(cowplot)
# trait <- "mQTL"
trait <- "haqtl"
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
annotations <- c("Repressed_Hoffman.to_hg38",
                 "Intron_UCSC.to_hg38",
                 "SuperEnhancer_Hnisz.to_hg38",
                 "UTR_3_UCSC.to_hg38",
                 "CTCF_Hoffman.to_hg38",
                 "DHS_peaks_Trynka.to_hg38",
                 "Enhancer_Andersson.to_hg38",
                 "PromoterFlanking_Hoffman.to_hg38",
                 "WeakEnhancer_Hoffman.to_hg38",
                 "FetalDHS_Trynka.to_hg38",
                 "Enhancer_Hoffman.to_hg38",
                 "DNaseI.to_hg38",
                 "Coding_UCSC.to_hg38",
                 "Promoter_UCSC.to_hg38",
                 "E081-DNase.macs2.narrowPeak.to_hg38",
                 "TSS_Hoffman.to_hg38",
                 "UTR_5_UCSC.to_hg38")
rows <- which(is.element(dat1$Annotation,annotations))
plot(dat1[rows,"Enrichment"],dat2[rows,"Enrichment"],pch = 20)
abline(a = 0,b = 1,col = "magenta",lty = "dotted")
dat1 <- dat1[rows,]
dat2 <- dat2[rows,]
# rows <- order(pmax(dat1$Enrichment,dat2$Enrichment))
# dat1 <- dat1[rows,]
# dat2 <- dat2[rows,]
# annotations <- dat1$Annotation
pdat <- rbind(data.frame(method = "SuSiE-topPC",
                         annotation = dat1$Annotation,
                         enrichment = dat1$Enrichment,
                         se         = dat1$Enrichment_SE,
                         stringsAsFactors = FALSE),
              data.frame(method = "fSuSiE",
                         annotation = dat2$Annotation,
                         enrichment = dat2$Enrichment,
                         se         = dat2$Enrichment_SE,
                         stringsAsFactors = FALSE))
pdat <- transform(pdat,
                  method     = factor(method),
                  annotation = factor(annotation,annotations),
                  enrichment = log2(enrichment))
print(ggplot(pdat,aes(x = enrichment,y = annotation,color = method,
                      xmin = enrichment - 2*se,
                      xmax = enrichment + 2*se)) +
      geom_point(shape = 20,size = 3) +
      geom_errorbarh(height = 0) +
      scale_color_manual(values = c("magenta","dodgerblue")) +
      scale_x_continuous(limits = c(-1.5,5),breaks = seq(-2,5)) +
      labs(x = "log2 enrichment",y = "",title = trait) +
      theme_cowplot(font_size = 10))
