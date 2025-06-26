rm(list=ls(
  
))
library(ggplot2)
library(tidyr)
library(AnnotationHub)
library(org.Hs.eg.db)
library(GenomicRanges)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(fsusieR)
library(dplyr)
library(data.table)
path= getwd()

cex=1
path= getwd()

obj_plot=list()
AD_GWAS <- fread("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/GWAS_sumstat.chr20_55439357_57610823.AD_Jansen_2021.tsv") 
obj_plot$AD_GWAS =AD_GWAS
AD_GWAS=obj_plot$AD_GWAS
qTLdata= fread(paste0("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/eQTL.chr20_55439357_57610823.tsv.tsv"))
obj_plot$qTLdata=qTLdata
qTLdata=obj_plot$qTLdata


chr=20
view_win=c(min(AD_GWAS$pos),max(AD_GWAS$pos))
pdat1= AD_GWAS 
dim(pdat1)
#pdat1$X.log10.P. 

pdat1$lg10p=-pnorm(-abs(pdat1$z), log.p = TRUE) / log(10)

t1= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start =pdat1$pos ,
                                                                 end = pdat1$pos )),
                data = matrix( (pdat1$lg10p) , nrow=1), genome = "hg19",
                ylim =c( min( (pdat1$lg10p)), max( (pdat1$lg10p))+0.2),
                type = "p", col = "black",  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
                background.title = "white",name="AD") ) # Change title color to black

plotTracks(t1)

pdat1CS =pdat1[which ( pdat1$CS1),]  
t2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = pdat1CS$pos   , 
                                                                 end = pdat1CS$pos   )),
                data = matrix( (pdat1CS$lg10p) , nrow=1), genome = "hg19", 
                ylim =c( min( (pdat1$lg10p)), max( (pdat1$lg10p))+0.2),
                type = "p", col = "black", cex=1.5,
                fill=  "royalblue",
                pch=c(24 ),# Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
                background.title = "white",name="AD") ) # Change title color to black
plotTracks(t2)
otAD <- OverlayTrack(trackList=list(    t1,  t2  ),
                     background.title = "white")

plotTracks( otAD , from= view_win[1],
            to= view_win[2])









## CASS4 panel ----- 



pdat2 =  qTLdata[which(qTLdata$study=="Mic_DeJager_eQTL"),]
#DLPFC also is consistent
#pdat2 = qTLdata[which(qTLdata$study=="DLPFC_DeJager_eQTL"),]

z_to_neglog10p <- function(z) {
  z <- abs(z)
  log_p <- pnorm(z, lower.tail = FALSE, log.p = TRUE)  # log(1 - Φ(z))
  log10_p <- log_p / log(10)
  -log10(2) - log10_p  # two-sided
} 
pdat2$X.log10.P. <- z_to_neglog10p(pdat2$z)


t1= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start =pdat2$pos ,
                                                                 end = pdat2$pos )),
                data = matrix(pdat2$X.log10.P. , nrow=1), genome = "hg19",
                ylim =c( min(pdat2$X.log10.P.), max(pdat2$X.log10.P.)+0.2),
                type = "p", col = "black",  # Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
                background.title = "white",name="CASS4") ) # Change title color to black

plotTracks(t1, from= view_win[1],
           to= view_win[2])

pdat2CS = pdat2[which(! is.na( pdat2$cs_coverage_0.95)),]
t2= ( DataTrack(range = GRanges(seqnames = chr, ranges = IRanges(start = pdat2CS$pos   , 
                                                                 end = pdat2CS$pos   )),
                data = matrix(pdat2CS$X.log10.P.  , nrow=1), genome = "hg19", 
                ylim =c( min(pdat2$X.log10.P.), max(pdat2$X.log10.P.)+0.2),
                type = "p", col = "black", cex=1.5,
                fill=  "royalblue",
                pch=c(24 ),# Use color column from df_plot
                track.margin = 0.05, # Reduce margin between track and title
                cex.title = 0.6,     # Reduce title size
                cex.axis = 0.6,      # Reduce axis text size
                col.axis = "black",  # Change axis color to black
                col.title = "black",cex.title = cex,
                rotation.title = 90,
                background.title = "white",name="AD") ) # Change title color to black


oteqTL <- OverlayTrack(trackList=list(    t1,  t2  ),
                       background.title = "white")

plotTracks( oteqTL, from= view_win[1],
            to= view_win[2])

plotTracks(list(otAD, oteqTL),
           
           from = 56407019, to =56433488)




# pip plot  -----

nrow(geno_raw)  -length( which ( geno_raw[,snp] == geno_raw[,best])) - length( which(is.na(geno_raw[,snp])))








#### meqtl -----



res <- readRDS("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/ROSMAP_mQTL.chr20_53859688_57519449.fsusie_mixture_normal_top_pc_weights.rds")
fsusie_obj_me = res$`chr20:53859688-57519449`$ROSMAP_DLPFC_mQTL$fsusie_result
rm(res) 
obj_plot$fsusie_obj_me = fsusie_obj_me 
res_me <- readRDS("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/ROSMAP_mQTL.chr20_53859688_57519449.fsusie_mixture_normal_top_pc_weights.input_data.rds")## to work from here
snp_names=attr( fsusie_obj_me$pip, "names")
pos_SNP_me <- as.numeric(sub("chr[0-9XY]+:([0-9]+):.*", "\\1", snp_names))






tt=  as.matrix(res_me$residual_X[[1]])
temp =cor( tt[,fsusie_obj_me$cs[[18]] ],
           tt[,-fsusie_obj_me$cs[[18]] ],
)
which.max(temp)
tt2= tt[,-fsusie_obj_me$cs[[18]] ]
plot(tt[,fsusie_obj_me$cs[[18]] ],
     tt2[,which.max(temp)],
     xlab="SNP in CS",
     ylab="Closest SNP to the SNP in CS")



 