# sinteractive -p mstephens -c 4 --mem=16G --time=10:00:00
# module load R/3.5.1
datadir <- "/project2/mstephens/fungen_xqtl/fsusie/mqtl/rosmap"
rds_files <- Sys.glob(paste0(datadir,"/ROSMAP_mQTL.chr*.fsusie_mixture_",
                             "normal_weights_db.rds"))
i <- 100
dat <- readRDS(rds_files[i])
dat <- dat[[1]]$ROSMAP_DLPFC_mQTL
# 
# > dat$region_info
# $grange
#   chrom    start      end
# 1    11 64525075 68955803
# 
# $region_name
# [1] "TADB_906"
#
# > head(dat$susie <- on <- top <- pc$variant <- names,n=3)
# [1] "chr11:64525133:G:A" "chr11:64525398:C:T" "chr11:64525649:T:G"
#
# dat$fsusie_result
#   $L_max
#   $L
#
# > dat$fsusie_result$N
# [1] 636
# > dat$fsusie_result$P
# [1] 17162
# 
# > dim(dat$fsusie_result$fitted_wc[[1]])
# [1] 17162  1024
#
