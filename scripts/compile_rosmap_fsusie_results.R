# TO DO: Explain here what this script is for, and how to use it.
#
# sinteractive -c 4 --mem=24G --time=120:00:00
# module load R/3.6.1
# R
# > .libPaths()[1]
# [1] "/home/pcarbo/R_libs_3_6"
options(stringsAsFactors = FALSE)

# TO DO: Add notes about these choices.
molecular_trait <- "mQTL"
analysis <- "ROSMAP_DLPFC_mQTL"
outfile <- "fsusie_ROSMAP_DLPFC_mQTL.RData"
datadir <- file.path("/project2/mstephens/fungen_xqtl/ftp_fgc_xqtl",
                     "analysis_result/finemapping_twas/fsusie")
fsusie_files <-
  Sys.glob(file.path(datadir,
    paste0("ROSMAP_",molecular_trait,
           ".*.fsusie_mixture_normal_top_pc_weights.rds")))

# Set up the data structures for storing the compiled results.
n <- length(fsusie_files)

# Repeat for each of the files to process.
cat("Compiling data from",n,"files:\n")
for (i in 1:n) {
  cat(i,"")
  dat <- readRDS(fsusie_files[i])
  dat <- dat[[1]][[analysis]]
}
cat("\n")

stop()

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
