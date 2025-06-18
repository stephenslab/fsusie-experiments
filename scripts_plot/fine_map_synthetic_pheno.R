rm(list=ls())
library(susieR)
#load data
res <- readRDS("/project2/mstephens/fungen_xqtl/CD2AP/data/ROSMAP_mQTL.chr6_44880827_48309905.fsusie_mixture_normal_top_pc_weights.rds")

fsusie_obj_me = res$`chr6:44880827-48309905`$ROSMAP_DLPFC_mQTL$fsusie_result
rm(res)
res_me <- readRDS("/project2/mstephens/fungen_xqtl/CD2AP/data/ROSMAP_mQTL.chr6_44880827_48309905.fsusie_mixture_normal_top_pc_weights.input_data_1.rds")
## to work from here



snp_names=attr( fsusie_obj_me$pip, "names")
pos_SNP_me <- as.numeric(sub("chr[0-9XY]+:([0-9]+):.*", "\\1", snp_names))

#this are the acutal data that were put into fsusie
Y= as.data.frame(res_me$residual_Y)
X=as.data.frame(res_me$residual_X)
pos = as.data.frame(res_me$Y_coordinates)#positions of the CpGs
# this the causal SNP displayed in figure 6A
fsusie_obj_me$cs[[13]]

#generate the synthetic pheno
#and other statistics
temp_pheno=list()
temp_pv=list()
temp_z=list()
temp_pheno_z=list()
for ( k in 1:ncol(Y)){

  fit_lm            = lm(Y[,k]~X[,fsusie_obj_me$cs[[13]]])
  temp_pheno[[k]]   = (summary(fit_lm)$coefficients[2,1] )*Y[,k]
  temp_pheno_z[[k]] = (summary(fit_lm)$coefficients[2,1]/summary(fit_lm)$coefficients[2,2] )*Y[,k]
  temp_pv[[k]]      = (summary(fit_lm)$coefficients[2,4] )
  temp_z[[k]]       = (summary(fit_lm)$coefficients[2,3] )
}


pv      = do.call(c, temp_pv)
zscore  = do.call(c, temp_z)
pheno   = Reduce("+",temp_pheno)
pheno_z = Reduce("+",temp_pheno_z)




m1= susie(y=pheno , X=as.matrix(X))
m2= susie(y=pheno_z , X=as.matrix(X))
m1$sets
m2$sets
# the fsusie result is
fsusie_obj_me$cs[[13]]

# correlation between first CS of susie and the CS of interest of fSuSiE
cor( X[ c(m2$sets$cs[[1]] , fsusie_obj_me$cs[[13]])])
#these snps seems to differs in only two individuals
plot(X[,m2$sets], X[, fsusie_obj_me$cs[[13]]],
     xlim= "FSuSiE SNP",
     ylim="SuSiE SNP")


#roughly transforming back the SNP data to the original format

SNP_f=   X[, fsusie_obj_me$cs[[13]]]

hist( SNP_f)
SNP_f[which(SNP_f<0) ]=0
SNP_f[which(SNP_f> 0& SNP_f<1) ]=1
SNP_f[which(SNP_f>1) ]=2
table(SNP_f)

#here the Manhattan plot for
#the fine map SNP cs each CpG
plot(pos[,2],-log10(pv),pch=21,
     xlab="genomic position")

#The region that Gao and his team fine mapped was actually quite large
# so in Fig 6A we focus on
plot(pos[,2],-log10(pv),pch=21,
     xlim=c(47400000,47700000),
     xlab="genomic position")




idx_zoom=which(pos[,2]> 47400000 & pos[,2] <47700000)


SNP_f
Y[,idx_zoom]
-log10(pv)[idx_zoom]

plot(zscore)



Y_zoom =Y[,idx_zoom]



library(ggplot2)
library(reshape2)
library(dplyr)


# showing cpG distrubition conditional on SNP value
# Step 0: Strip off any non-feature column (e.g., Sample)
Y_features_only <- Y_zoom[, sapply(Y_zoom, is.numeric)]  # keep only numeric columns
feature_names <- colnames(Y_features_only)
genomic_positions <- pos[idx_zoom, 2]
stopifnot(length(feature_names) == length(genomic_positions))
Y_features_only$Sample <- rownames(Y_zoom)
Y_long <- melt(Y_features_only, id.vars = "Sample", variable.name = "Feature", value.name = "Value")
n <- nrow(Y_zoom)
p <- length(feature_names)
Y_long$SNP_f <- factor(rep(SNP_f, times = p))
position_df <- data.frame(Feature = feature_names, GenomicPos = genomic_positions)
Y_long <- left_join(Y_long, position_df, by = "Feature")



ggplot(Y_long, aes(x = GenomicPos, y = Value, fill = SNP_f, group = interaction(GenomicPos, SNP_f))) +
  geom_boxplot(  ) +
  labs(x = "Genomic Position", y = "Expression / Value", fill = "SNP Genotype") +
  theme_minimal() +
  ylim(c(-0.1,0.1))+
  theme(legend.position = "top")

# The problem is that these CpG are so far apart that it is kind of hard to properly distinguish


#Putting the CpG boxplot next to each other
ggplot(Y_long, aes(x = Feature, y = Value, fill = SNP_f)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.3, lwd = 0.3) +
  labs(x = "Feature", y = "Expression / Value", fill = "SNP Genotype") +
  theme_minimal() +
  ylim(c(-0.1,0.1))+

  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "top")

