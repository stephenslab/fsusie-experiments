rm(list = ls())

library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

## =========================
## PATHS
## =========================
result_path <- "/home/wdenault/fsusi_simu/correlated_sim/permutation_scATAC_results/"
lf <- list.files(result_path,
                 pattern = "^result_scATAC_.*\\.RData$",
                 full.names = TRUE)

stopifnot(length(lf) > 0)

## =========================
## HELPERS
## =========================

get_cs_sizes <- function(cs) {
  if (is.null(cs) || length(cs) == 0)
    return(integer(0))
  unlist(lapply(cs, length), use.names = FALSE)
}

safe_load <- function(f) {
  e <- new.env()
  load(f, envir = e)
  e$out
}

## =========================
## LOAD ALL RESULTS
## =========================
all_res <- map(lf, safe_load)

## =========================
## FLATTEN
## =========================
df <- map_dfr(all_res, function(x) {

  param <- as.list(x$param)

  map_dfr(x$results, function(r) {

    tibble(
      nullweight      = as.numeric(param$nullweight),
      thresh_lowcount = as.numeric(param$thresh_lowcount),
      region          = r$region,

      mfsusie_cs            = list(get_cs_sizes(r$mfsusie)),
      mfsusie_corsmall_cs   = list(get_cs_sizes(r$mfsusie_corsmall)),
      fsusie_cs             = list(get_cs_sizes(r$fsusie)),
      fsusie_corsmall_cs    = list(get_cs_sizes(r$fsusie_corsmall))
    )
  })
})

## =========================
## LONG FORMAT
## =========================
df_long <- df %>%
  pivot_longer(
    cols = ends_with("_cs"),
    names_to = "method",
    values_to = "cs_sizes"
  ) %>%
  mutate(
    method = recode(method,
                    mfsusie_cs          = "mfsusie",
                    mfsusie_corsmall_cs = "mfsusie_corsmall",
                    fsusie_cs           = "fsusie",
                    fsusie_corsmall_cs  = "fsusie_corsmall"
    ),
    n_cs   = map_int(cs_sizes, length),
    has_cs = n_cs > 0
  ) %>%
  unnest(cs_sizes, keep_empty = TRUE) %>%
  rename(cs_size = cs_sizes)

## =========================
## SUMMARY TABLE
## =========================
summary_tbl <- df_long %>%
  group_by(method, nullweight, thresh_lowcount) %>%
  summarise(
    n_runs = n_distinct(paste(region, row_number())),
    prop_nonempty = mean(has_cs),
    mean_cs_size  = mean(cs_size, na.rm = TRUE),
    median_cs_size = median(cs_size, na.rm = TRUE),
    max_cs_size   = max(cs_size, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_tbl)

## =========================
## NULL CALIBRATION PLOTS
## =========================

## CS size distribution
ggplot(df_long %>% filter(!is.na(cs_size)),
       aes(x = cs_size, fill = method)) +
  geom_histogram(bins = 40, alpha = 0.6, position = "identity") +
  facet_grid(nullweight ~ thresh_lowcount) +
  scale_x_continuous(trans = "log1p") +
  theme_bw() +
  scale_x_log10()+
  labs(
    title = "Credible set sizes under permutation null (scATAC)",
    x = "CS size (log1p)",
    y = "Count"
  )

## Probability of any CS
ggplot(df_long %>%
         distinct(method, nullweight, thresh_lowcount, region, has_cs),
       aes(x = factor(thresh_lowcount),
           y = has_cs,
           color = factor(nullweight),
           group = nullweight)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun = mean, geom = "point") +
  facet_wrap(~ method, ncol = 2) +
  theme_bw() +
  labs(
    y = "P(any CS | null)",
    x = "thresh_lowcount",
    color = "nullweight"
  )
