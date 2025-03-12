path="D:/Document/Serieux/Travail/Data_analysis_and_papers/GTEX_analysis_Fsusie"
source(paste0(path,"/code/plot_all_effect_log.R"))
source(paste0(path,"/code/plot_log.R"))

load(paste0(path,"/output/local/ENSG00000080824.csv.gz.RData"))
plot_fsusie_log (out = out, log_effect = FALSE   ) 


library(fsusieR)

plot_susiF_pip(out$res, pos_SNP = as.numeric(out$info_SNP$POS))+
  
  geom_rect(aes(xmin = out$locus[1], 
                xmax =out$locus[2],
                ymin = -0.01,
                ymax = 1.02), 
            alpha = 0.0, color = "red")


#load(paste0(path,"/output/df_brain_cortex_junction.RData")) 
fsusie_log_plot(out$res,chr = paste0("chr",out$chr),
                pos0 = out$locus[1],pos1 = out$locus[2],
                out$X,out$Y,snp_info = out$info_SNP,cs = 1,
                effect_log=TRUE,
                log1p_count=TRUE )

fsusie_log_plot(out$res,chr = paste0("chr",out$chr),
                pos0 = out$locus[1],pos1 = out$locus[2],
                out$X,out$Y,snp_info = out$info_SNP,cs = 2,
                effect_log=TRUE,
                log1p_count=TRUE )
 
