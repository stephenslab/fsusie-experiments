# TO DO: Explain here what this script is for, and how to use it.

# module load R/3.6.1
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
                            snps         = rep(0,n),
                            stringsAsFactors = FALSE))

# Repeat for each of the files to process.
# for (i in 1:n) {
for (i in 3372) {
  cat(susie_files[i],"\n")
  dat <- readRDS(susie_files[i])
  dat <- dat[[1]][[analysis]]
  m   <- length(dat$pip)

  # Get the region info.
  susie_rnaseq$regions[i,"region_name"]  <- dat$region_info$region_name
  susie_rnaseq$regions[i,"chr"]          <- dat$region_info$region_coord$chrom
  susie_rnaseq$regions[i,"coord_start"]  <- dat$region_info$region_coord$start
  susie_rnaseq$regions[i,"coord_end"]    <- dat$region_info$region_coord$end
  susie_rnaseq$regions[i,"grange_chr"]   <- dat$region_info$grange$chrom
  susie_rnaseq$regions[i,"grange_start"] <- dat$region_info$grange$start
  susie_rnaseq$regions[i,"grange_end"]   <- dat$region_info$grange$end
  susie_rnaseq$regions[i,"snps"]         <- m

  # Get the PIPs.
  susie_rnaseq$pips[[i]] <-
    data.frame(region = dat$region_info$region_name,
               id     = names(dat$pip),
               pos    = sapply(strsplit(names(dat$pip),":"),"[",2),
               pip    = dat$pip,
               stringsAsFactors = FALSE)
  rownames(susie_rnaseq$pips[[i]]) <- NULL
  
  stop()
}

# susie_get_cs with coverage = 0.95
# min_abs_corr > 0.5 OR median_abs_corr > 0.8

# > dat$region_info$region_name
# [1] "ENSG00000105369" "ENSG00000105369"
# > dat$region_info$region_coord
#   chrom    start      end
# 1    19 41877278 41877279
# > dat$region_info$grange
#   chrom    start      end
# 1    19 35640000 45640000
# > head(dat$pip,n = 3)
#    chr19:35641273:C:T chr19:35641582:ACTC:A    chr19:35642020:T:C
#             0.0003640             0.0002389             0.0002920

