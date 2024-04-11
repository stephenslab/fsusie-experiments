# TO DO: Explain here what this script is for, and how to use it.

# sinteractive -c 4 --mem=16G --time=24:00:00
# module load R/3.6.1
# R
# > .libPaths()[1]
# [1] "/home/pcarbo/R_libs_3_6"

# This is one of the QTL results that Gao suggested to prioritize.
# (Does Hao agree?)
analysis <- "AC_DeJager_eQTL"

datadir <- file.path("/project2/mstephens/fungen_xqtl/ftp_fgc_xqtl",
                     "analysis_result/finemapping_twas",
                     "susie_twas_export_new")
susie_files <-
  Sys.glob(file.path(datadir,"Fungen_xQTL.ENSG*.cis_results_db.export.rds"))

# Set up the data structures for storing the compiled results.
n <- length(susie_files)
susie_rnaseq <-
  list(pips = vector("list",n),
       cs   = vector("list",n),
       regions = data.frame(region_name  = rep("",n),
                            chr          = rep(0,n),
                            coord_start  = rep(0,n),
                            coord_end    = rep(0,n),
                            grange_chr   = rep(0,n),
                            grange_start = rep(0,n),
                            grange_end   = rep(0,n),
                            num_snps     = rep(0,n),
                            num_cs       = rep(0,n),
                            stringsAsFactors = FALSE))

# Repeat for each of the files to process.
for (i in 1:n) {
  cat(i,"")
  dat <- readRDS(susie_files[i])
  dat <- dat[[1]][[analysis]]
  if (is.null(dat$pip))
    next

  # Get the number of SNPs.
  m <- length(dat$pip)
  
  # Get the number of credible sets (CSs).
  if (nrow(dat$top_loci) > 0)
    num_cs <- max(dat$top_loci$cs_coverage_0.95)
  else
    num_cs <- 0
  
  # Get the region info.
  susie_rnaseq$regions[i,"region_name"]  <- dat$region_info$region_name
  susie_rnaseq$regions[i,"chr"]          <- dat$region_info$region_coord$chrom
  susie_rnaseq$regions[i,"coord_start"]  <- dat$region_info$region_coord$start
  susie_rnaseq$regions[i,"coord_end"]    <- dat$region_info$region_coord$end
  susie_rnaseq$regions[i,"grange_chr"]   <- dat$region_info$grange$chrom
  susie_rnaseq$regions[i,"grange_start"] <- dat$region_info$grange$start
  susie_rnaseq$regions[i,"grange_end"]   <- dat$region_info$grange$end
  susie_rnaseq$regions[i,"num_snps"]     <- m
  susie_rnaseq$regions[i,"num_cs"]       <- num_cs

  # Get the PIPs.
  susie_rnaseq$pips[[i]] <-
    data.frame(region = dat$region_info$region_name,
               id     = names(dat$pip),
               pos    = sapply(strsplit(names(dat$pip),":"),"[",2),
               pip    = dat$pip,
               stringsAsFactors = FALSE)
  rownames(susie_rnaseq$pips[[i]]) <- NULL

  # Get the credible sets (CSs). The CSs are obtain using susie_get_cs()
  # with coverage = 0.95, then filtering out CSs with min_abs_corr <
  # 0.5 and median_abs_corr < 0.8.
  if (num_cs == 0)
    cs <- as.character(NA)
  else {
    cs <- data.frame(region    = dat$region_info$region_name,
                     id        = dat$top_loci$variant_id,
                     betahat   = dat$top_loci$betahat,
                     sebetahat = dat$top_loci$sebetahat,
                     maf       = dat$top_loci$maf,
                     pip       = dat$top_loci$pip,
                     cs        = dat$top_loci$cs_coverage_0.95)
    cs$cs[cs$cs == 0] <- NA
  }
  susie_rnaseq$cs[[i]] <- cs
}
