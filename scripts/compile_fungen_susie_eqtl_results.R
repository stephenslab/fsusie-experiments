# TO DO: Explain here what this script is for, and how to use it.
#
# sinteractive -c 4 --mem=16G --time=120:00:00
# module load R/3.6.1
# R
# > .libPaths()[1]
# [1] "/home/pcarbo/R_libs_3_6"
options(stringsAsFactors = FALSE)

# This is one of the QTL results that Gao suggested to prioritize.
# (Does Hao agree?)
#
# Also I want this to match up these results with the fsusie results
# on methylation ("ROSMAP_DLPFC_mQTL").
analysis <- "Inh_mega_eQTL"
outfile <- "susie_Inh_mega_eQTL.RData"
datadir <- file.path("/project2/mstephens/fungen_xqtl/ftp_fgc_xqtl",
                     "analysis_result/finemapping_twas",
                     "susie_twas_export_new")
susie_files <-
  Sys.glob(file.path(datadir,"Fungen_xQTL.ENSG*.cis_results_db.export.rds"))

# Set up the data structures for storing the compiled results.
n       <- length(susie_files)
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
                      num_cs       = rep(0,n))

# Repeat for each of the files to process.
cat("Compiling data from",n,"files:\n")
for (i in 1:n) {
  cat(i,"")
  dat <- readRDS(susie_files[i])
  dat <- dat[[1]][[analysis]]

  # If there are no SNPs, skip the region.
  if (is.null(dat$pip))
    next

  # Get the number of SNPs.
  m <- length(dat$pip)
  
  # Get the number of credible sets (CSs).
  if (is.null(dat$top_loci))
    num_cs <- 0
  else if (nrow(dat$top_loci) == 0)
    num_cs <- 0
  else
    num_cs <- max(dat$top_loci$cs_coverage_0.95)
  
  # Get the region info.
  regions[i,"region_name"]  <- dat$region_info$region_name[1]
  regions[i,"chr"]          <- dat$region_info$region_coord$chrom
  regions[i,"coord_start"]  <- dat$region_info$region_coord$start
  regions[i,"coord_end"]    <- dat$region_info$region_coord$end
  regions[i,"grange_chr"]   <- dat$region_info$grange$chrom
  regions[i,"grange_start"] <- dat$region_info$grange$start
  regions[i,"grange_end"]   <- dat$region_info$grange$end
  regions[i,"num_snps"]     <- m
  regions[i,"num_cs"]       <- num_cs

  # Get the PIPs.
  pips[[i]] <-
    data.frame(region = dat$region_info$region_name[1],
               id     = names(dat$pip),
               pos    = sapply(strsplit(names(dat$pip),":"),"[",2),
               pip    = dat$pip)
  rownames(pips[[i]]) <- NULL

  # Get the credible sets (CSs). The CSs are obtain using susie_get_cs()
  # with coverage = 0.95, then filtering out CSs with min_abs_corr <
  # 0.5 and median_abs_corr < 0.8.
  if (num_cs == 0)
    res <- as.character(NA)
  else {
    res <- data.frame(region    = dat$region_info$region_name[1],
                          id        = dat$top_loci$variant_id,
                          betahat   = dat$top_loci$betahat,
                          sebetahat = dat$top_loci$sebetahat,
                          maf       = dat$top_loci$maf,
                          pip       = dat$top_loci$pip,
                          cs        = dat$top_loci$cs_coverage_0.95)
    res$cs[res$cs == 0] <- NA
  }
  cs[[i]] <- res
}
cat("\n")
regions <- transform(regions,
                     chr        = as.numeric(chr),
                     grange_chr = as.numeric(grange_chr))
rm(datadir,susie_files)
rm(n,m,i,dat,res,num_cs)

# Save the final data structure to an RDS file.
save(file = outfile,list = c("regions","pips","cs"))

stop()

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

stop()

# Combine the CS results into a single data frame.
i  <- which(!is.na(cs))
cs <- cs[i]
cs <- do.call(rbind,cs)
i  <- which(!is.na(cs$cs))
cs <- cs[i,]
# cs <- transform(cs)

# Combine the PIPs into a single data frame.
pips <- do.call(rbind,pips)

# Save the final data structure to an RDS file.
save(file = outfile,list = c("regions","pips","cs"))
