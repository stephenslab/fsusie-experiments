
rm(list = ls())
set.seed(100030)

library(fsusieR)
library(susieR)

input_path  <- "/project2/mstephens/wdenault/GTEX_analysis_Fsusie/result_Gao/"
result_file <- "/home/wdenault/fsusi_simu/correlated_sim/results//result_perm_30.RData"

lf <- list.files(input_path)

res_job <- vector("list", 20)

for (k in seq_len(20)) {

  id <- sample(seq_along(lf), 1)
  load(paste0(input_path, lf[id]))

  X_perm <- as.matrix(out$X[sample(nrow(out$X)), ])

  tt <- fsusieR::susiF(
    Y = log1p(out$Y),
    X = X_perm,
    L = 2,
    nullweight = 0.100,
    post_processing = "none",
    quantile_trans = TRUE,
    thresh_lowcount = 0.500
  )

  tt_small <- fsusieR::susiF(
    Y = log1p(out$Y),
    X = X_perm,
    L = 2,
    nullweight = 0.100,
    post_processing = "none",
    quantile_trans = TRUE,
    cor_small = TRUE,
    thresh_lowcount = 0.500
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
  param = c(nullweight = 0.100, thresh_lowcount = 0.500),
  job_id = 30
)

save(out_res, file = result_file)
}



