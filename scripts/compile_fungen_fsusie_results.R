# This script compiles the fsusie results into data structures for
# more convenient exploratory analyses.
#
# sinteractive -c 4 --mem=96G --time=120:00:00
# module load R/3.6.1
# R
# > .libPaths()[1]
# [1] "/home/pcarbo/R_libs_3_6"
library(tools)
options(stringsAsFactors = FALSE)
#
# molecular_trait <- "mQTL"
# analysis <- "ROSMAP_DLPFC_mQTL"
# outfile <- "fsusie_ROSMAP_DLPFC_mQTL.RData"
#
molecular_trait <- "haQTL"
analysis <- "ROSMAP_DLPFC_haQTL"
outfile <- "fsusie_ROSMAP_DLPFC_haQTL.RData"
datadir <- file.path("/project2/mstephens/fungen_xqtl/ftp_fgc_xqtl",
                     "analysis_result/finemapping_twas/fsusie")
fsusie_files <-
  Sys.glob(file.path(datadir,
    paste0("ROSMAP_",molecular_trait,
           ".*.fsusie_mixture_normal_top_pc_weights.rds")))

# Set up the data structures for storing the compiled results.
n       <- length(fsusie_files)
pips    <- vector("list",n)
cs      <- vector("list",n)
regions <- data.frame(region_name  = rep("",n),
                      chr          = rep(0,n),
                      coord_start  = rep(0,n),
                      coord_end    = rep(0,n),
                      grange_chr   = rep(0,n),
                      grange_start = rep(0,n),
                      grange_end   = rep(0,n),
                      num_snps     = rep(0,n),
                      num_cs       = rep(0,n),
                      L            = rep(0,n),
                      Lmax         = rep(0,n))

# Repeat for each of the files to process.
cat("Compiling data from",n,"files:\n")
for (i in 1:n) {
  cat(i,"")
  tryCatch
  dat <- tryCatch(readRDS(fsusie_files[i]),
                  error = function (e) NULL)
  if (is.null(dat))
    next
  dat <- dat[[1]][[analysis]]

  # Get the number of SNPs.
  # If there are no SNPs, skip the region.
  num_snps <- length(dat$fsusie_summary$variant_names)
  if (num_snps == 0)
    next
  
  # Get the number of credible sets (CSs).
  top_loci <- dat$fsusie_summary$top_loci
  if (is.null(top_loci))
    num_cs <- 0
  else if (nrow(top_loci) == 0)
    num_cs <- 0
  else
    num_cs <- max(top_loci$cs_coverage_0.95)
  
  # Get the region info.
  regions[i,"region_name"]  <- dat$region_info$region_name
  regions[i,"chr"]          <- dat$region_info$region_coord$chrom
  regions[i,"coord_start"]  <- dat$region_info$region_coord$start
  regions[i,"coord_end"]    <- dat$region_info$region_coord$end
  regions[i,"grange_chr"]   <- dat$region_info$grange$chrom
  regions[i,"grange_start"] <- dat$region_info$grange$start
  regions[i,"grange_end"]   <- dat$region_info$grange$end
  regions[i,"num_snps"]     <- num_snps
  regions[i,"num_cs"]       <- num_cs
  regions[i,"L"]            <- dat$fsusie_result$L
  regions[i,"Lmax"]         <- dat$fsusie_result$L_max
  
  # Get the PIPs.
  pips[[i]] <-
    data.frame(region = dat$region_info$region_name[1],
               id = names(dat$fsusie_result$pip),
               pos = sapply(strsplit(names(dat$fsusie_result$pip),":"),"[",2),
               pip = as.numeric(dat$fsusie_result$pip))
  rownames(pips[[i]]) <- NULL
  
  # Get the credible sets (CSs). The CSs are obtain using susie_get_cs()
  # with coverage = 0.95, then filtering out CSs with min_abs_corr <
  # 0.5 and median_abs_corr < 0.8.
  if (num_cs == 0)
    res <- as.character(NA)
  else {
    res <- data.frame(region    = dat$region_info$region_name[1],
                      id        = top_loci$variant_id,
                      maf       = top_loci$maf,
                      pip       = top_loci$pip,
                      cs        = top_loci$cs_coverage_0.95)
    res$cs[res$cs == 0] <- NA
  }
  cs[[i]] <- res
}
cat("\n")
regions <- transform(regions,
                     chr        = as.numeric(chr),
                     grange_chr = as.numeric(grange_chr))

# Consider using these results later:
#
#   dat$fsusie_result$fitted_wc
#     $fitted_wc
#     $cred_band
#     $est_pi
#

# Remove the regions that have no SNPs.
i       <- which(regions$num_snps > 0)
regions <- regions[i,]
pips    <- pips[i]
cs      <- cs[i]

# Sort the regions by chromosome and by position along the chromosome.
i       <- with(regions,order(chr,coord_start))
regions <- regions[i,]
pips    <- pips[i]
cs      <- cs[i]

# Combine the CS results into a single data frame.
i  <- which(!is.na(cs))
cs <- cs[i]
cs <- do.call(rbind,cs)
i  <- which(!is.na(cs$cs))
cs <- cs[i,]
cs <- transform(cs,
                region = factor(region),
                cs     = factor(cs))

# Add the region names to the PIPs data structure for easier lookup.
names(pips) <- regions$region_name

# Convert the SNP base-pair positions into numbers.
n <- length(pips)
for (i in 1:n) {
  pips[[i]] <- transform(pips[[i]],pos = as.numeric(pos))
}

# Save the final data structure to an RDS file.
save(file = outfile,list = c("regions","pips","cs"))
resaveRdaFiles(outfile)
