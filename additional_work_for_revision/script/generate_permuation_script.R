rm(list = ls())

## =========================
## USER SETTINGS
## =========================
input_path  <- "/project2/mstephens/wdenault/GTEX_analysis_Fsusie/result_Gao/"
output_path <- "/home/wdenault/fsusi_simu/script/additional_simulations/permutation_job/"
result_path <- "/home/wdenault/fsusi_simu/correlated_sim/results/"

nullweights <- c(0.1, 0.5)
thresh_lowcounts <- c(0.05, 0.1, 0.5)

n_rep_job <- 20      # replicates INSIDE each job
n_jobs_rep <- 10     # number of jobs per parameter combo

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(result_path, recursive = TRUE, showWarnings = FALSE)

lf <- list.files(input_path)

## =========================
## JOB GENERATION
## =========================
job_id <- 1

for (null_w in nullweights) {
  for (thres_c in thresh_lowcounts) {
    for (job_rep in seq_len(n_jobs_rep)) {

      job_file    <- paste0(output_path, "/job_", job_id, ".R")
      result_file <- paste0(result_path, "/result_perm_", job_id, ".RData")

      script <- sprintf(
        '
rm(list = ls())
set.seed(%d)

library(fsusieR)
library(susieR)

input_path  <- "%s"
result_file <- "%s"

lf <- list.files(input_path)

res_job <- vector("list", %d)

for (k in seq_len(%d)) {

  id <- sample(seq_along(lf), 1)
  load(paste0(input_path, lf[id]))

  X_perm <- as.matrix(out$X[sample(nrow(out$X)), ])

  tt <- fsusieR::susiF(
    Y = log1p(out$Y),
    X = X_perm,
    L = 2,
    nullweight = %.3f,
    post_processing = "none",
    quantile_trans = TRUE,
    thresh_lowcount = %.3f
  )

  tt_small <- fsusieR::susiF(
    Y = log1p(out$Y),
    X = X_perm,
    L = 2,
    nullweight = %.3f,
    post_processing = "none",
    quantile_trans = TRUE,
    cor_small = TRUE,
    thresh_lowcount = %.3f
  )

  lol <- susieR::susie(
    y = log1p(apply(out$Y / out$size_factor_local, 1, sum)),
    X = X_perm
  )

  res_job[[k]] <- list(
    fsusie = tt$cs,
    fsusie_corsmall = tt_small$cs,
    susie = lol$sets,
    file_used = lf[id]
  )
  out_res <- list(
  results = res_job,
  param = c(nullweight = %.3f, thresh_lowcount = %.3f),
  job_id = %d
)

save(out_res, file = result_file)
}


',
      ## RNG seed (unique per job)
      100000 + job_id,

      input_path,
      result_file,

      n_rep_job,
      n_rep_job,

      null_w, thres_c,
      null_w, thres_c,

      null_w, thres_c,

      job_id
      )

writeLines(script, con = job_file)
job_id <- job_id + 1
    }
  }
}

cat("Generated", job_id - 1, "job scripts in", output_path, "\n")
