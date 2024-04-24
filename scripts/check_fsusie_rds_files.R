# sinteractive -c 4 --mem=96G --time=120:00:00
# module load R/3.6.1
# R
datadir <- file.path("/project2/mstephens/fungen_xqtl/ftp_fgc_xqtl",
                     "analysis_result/finemapping_twas/fsusie")
molecular_trait <- "mQTL"
fsusie_files <-
  Sys.glob(file.path(datadir,
    paste0("ROSMAP_",molecular_trait,
           ".*.fsusie_mixture_normal_top_pc_weights.rds")))
n <- length(fsusie_files)
corrupted <- rep(FALSE,n)
cat("Checking",n,"files:\n")
for (i in 1:n) {
  cat("",i)
  dat <- tryCatch(readRDS(fsusie_files[i]),
                  error = function (e) NULL)
  if (is.null(dat)) {
    cat("*")
    corrupted[i] <- TRUE
  }
}
cat("\n")
