rm(list = ls())

## =========================
## USER SETTINGS
## =========================
input_path  <- "/project2/mstephens/wdenault/anjing_data/data/"
job_path    <- "/home/wdenault/fsusi_simu/script/additional_simulations/permutation_scATAC/"
result_path <- "/home/wdenault/fsusi_simu/correlated_sim/permutation_scATAC_results/"

nullweights <- c(0.1, 0.5)
thresh_lowcounts <- c(0.01, 0.05, 0.1, 0.5)

n_jobs_rep <- 10      # jobs per parameter combo
n_rep_job  <- 10      # permutations inside each job
filter_ncell <- 0

dir.create(job_path, recursive = TRUE, showWarnings = FALSE)
dir.create(result_path, recursive = TRUE, showWarnings = FALSE)

## =========================
## JOB GENERATION
## =========================
job_id <- 1

for (null_w in nullweights) {
  for (thres_c in thresh_lowcounts) {
    for (job_rep in seq_len(n_jobs_rep)) {

      job_file    <- paste0(job_path, "/job_", job_id, ".R")
      result_file <- paste0(result_path, "/result_scATAC_", job_id, ".RData")

      script <- sprintf(
        '
rm(list = ls())
set.seed(%d)

library(mvf.susie.alpha)
library(fsusieR)
library(susieR)

input_path  <- "%s"
result_file <- "%s"
filter_ncell <- %d

## -------------------------
## Load data
## -------------------------
ncells <- t(read.delim(paste0(input_path,
  "kellis_snatac.combined_genome_ncell.tss6.filtered_84samples.6celltype.txt"),
  row.names=1))
rownames(ncells) <- gsub("\\\\.", "-", rownames(ncells))

all_count <- t(read.delim(paste0(input_path,
  "updated.kellis_snatac.combined_genome_reads.tss6.filtered_84samples.6celltype.txt"),
  row.names=1))
rownames(all_count) <- gsub("\\\\.", "-", rownames(all_count))

counts <- readRDS(paste0(input_path,
  "xiong_atact_seq_multi_trait.169regions.51.2kb_expanded.binned_coverage_count.reordered_1mb_maf005.rds"))

## -------------------------
## Storage
## -------------------------
res_job <- vector("list", %d)

for (rep_id in seq_len(%d)) {

  reg_num <- sample(seq_along(counts), 1)
  tcount  <- counts[[reg_num]]
  X <- tcount$filtered_geno

  Y_f <- list()

  for (k in seq_len(length(tcount) - 1)) {

    tt <- tcount[[k]]$pheno
    common <- intersect(rownames(tt), rownames(ncells))

    nc_ord  <- ncells[common, , drop=FALSE]
    ac_ord  <- all_count[common, , drop=FALSE]

    avg_read <- ac_ord[,k] / nc_ord[,k]
    scaling  <- avg_read / mean(avg_read, na.rm=TRUE)

    temp_Y  <- tt * (1 / nc_ord[,k])
    Y_f[[k]] <- log1p(temp_Y / scaling)
  }

  ## NA handling
  for (i in seq_len(nrow(X))) {
    for (k in seq_along(Y_f)) {
      if (is.na(ncells[i,k])) Y_f[[k]][i,] <- NA
    }
  }

  if (filter_ncell > 0) {
    for (i in seq_len(nrow(X))) {
      for (k in seq_along(Y_f)) {
        if (ncells[i,k] < filter_ncell || is.na(ncells[i,k]))
          Y_f[[k]][i,] <- NA
      }
    }
  }

  ## Genotype cleanup
  for (k in seq_len(ncol(X))) {
    if (anyNA(X[,k]))
      X[is.na(X[,k]),k] <- median(X[,k], na.rm=TRUE)
  }

  maf <- apply(X,2,mean)/2
  X <- X[, maf >= 0.1 & maf <= 0.9, drop=FALSE]

  Y <- list(Y_u=NULL, Y_f=Y_f)

  X_perm <- as.matrix(X[sample(nrow(X)), ])

  pos <- replicate(6, 1:1024, simplify=FALSE)

  ## -------------------------
  ## Methods
  ## -------------------------
  tt_m <- multfsusie(
    Y=Y, X=X_perm, pos=pos, L=3,
    cor_small=FALSE,
    post_processing="none",
    quantile_trans=TRUE,
    thresh_lowcount=%.3f
  )

  tt_msmall <- multfsusie(
    Y=Y, X=X_perm, pos=pos, L=3,
    cor_small=TRUE,
    nullweight=%.3f,
    post_processing="none",
    quantile_trans=TRUE,
    thresh_lowcount=%.3f
  )

  id <- sample(seq_along(Y_f),1)

  tt <- fsusieR::susiF(
    Y=log1p(Y_f[[id]]), X=X_perm, L=2,
    nullweight=%.3f,
    post_processing="none",
    quantile_trans=TRUE,
    thresh_lowcount=%.3f
  )

  tt_small <- fsusieR::susiF(
    Y=log1p(Y_f[[id]]), X=X_perm, L=2,
    nullweight=%.3f,
    cor_small=TRUE,
    post_processing="none",
    quantile_trans=TRUE,
    thresh_lowcount=%.3f
  )

  res_job[[rep_id]] <- list(
    mfsusie = tt_m$cs,
    mfsusie_corsmall = tt_msmall$cs,
    fsusie = tt$cs,
    fsusie_corsmall = tt_small$cs,
    region = reg_num
  )

out <- list(
  results = res_job,
  param = c(nullweight=%.3f, thresh_lowcount=%.3f),
  job_id = %d
)

save(out, file=result_file)
}

',
      ## seed
      100000 + job_id,

      input_path,
      result_file,
      filter_ncell,

      n_rep_job,
      n_rep_job,

      thres_c,
      null_w, thres_c,
      null_w, thres_c,
      null_w, thres_c,

      null_w, thres_c,
      job_id
      )

writeLines(script, job_file)
job_id <- job_id + 1
    }
  }
}

cat("Generated", job_id - 1, "scATAC job scripts in", job_path, "\n")
