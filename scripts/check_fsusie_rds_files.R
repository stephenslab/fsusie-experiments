datadir <- file.path("/project2/mstephens/fungen_xqtl/ftp_fgc_xqtl",
                     "analysis_result/finemapping_twas/fsusie")
molecular_trait <- "mQTL"
fsusie_files <-
  Sys.glob(file.path(datadir,
    paste0("ROSMAP_",molecular_trait,
           ".*.fsusie_mixture_normal_top_pc_weights.rds")))
# *** TESTING ***
fsusie_files <- fsusie_files[1:100]
n <- length(fsusie_files)
corrupted <- rep(FALSE,n)
cat("Checking",n,"files:\n")
for (i in 1:n) {
  cat(i,"")
  dat <- tryCatch(readRDS(fsusie_files[i]),
                  error = function (e) 123)
  if (is.null(dat))
    corrupted[i] <- TRUE
}
