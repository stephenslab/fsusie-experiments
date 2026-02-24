rm(list = ls())

library(dplyr)
library(purrr)
library(tidyr)

## =========================
## PATHS
## =========================
result_path <- "/home/wdenault/fsusi_simu/correlated_sim/results/"
lf <- list.files(result_path, pattern = "^result_perm_.*\\.RData$", full.names = TRUE)

stopifnot(length(lf) > 0)

## =========================
## HELPER FUNCTIONS
## =========================

## Extract CS sizes from a susie-style CS object
get_cs_sizes <- function(cs) {
  if (is.null(cs) || length(cs) == 0) {
    return(integer(0))
  }
  unlist(lapply(cs, length), use.names = FALSE)
}

## Safe loader
load_result <- function(f) {
  e <- new.env()
  load(f, envir = e)
  e$out_res
}

## =========================
## LOAD & UNNEST RESULTS
## =========================

all_res <- map(lf, load_result)

df <- map_dfr(all_res, function(x) {

  param <- as.list(x$param)

  map_dfr(x$results, function(r) {

    tibble(
      nullweight       = as.numeric(param$nullweight),
      thresh_lowcount  = as.numeric(param$thresh_lowcount),

      fsusie_cs_sizes =
        list(get_cs_sizes(r$fsusie)),

      fsusie_corsmall_cs_sizes =
        list(get_cs_sizes(r$fsusie_corsmall))
    )
  })
})

## =========================
## LONG FORMAT
## =========================

df_long <- df %>%
  pivot_longer(
    cols = c(fsusie_cs_sizes, fsusie_corsmall_cs_sizes),
    names_to = "method",
    values_to = "cs_sizes"
  ) %>%
  mutate(
    method = recode(
      method,
      fsusie_cs_sizes = "fsusie",
      fsusie_corsmall_cs_sizes = "fsusie_corsmall"
    ),
    n_cs = map_int(cs_sizes, length),
    has_cs = n_cs > 0
  ) %>%
  unnest(cs_sizes, keep_empty = TRUE) %>%
  rename(cs_size = cs_sizes)

## =========================
## SUMMARY STATISTICS
## =========================

summary_tbl <- df_long %>%
  group_by(method, nullweight, thresh_lowcount) %>%
  summarise(
    n_runs = n_distinct(row_number()),
    prop_nonempty = mean(has_cs),
    mean_cs_size = mean(cs_size, na.rm = TRUE),
    median_cs_size = median(cs_size, na.rm = TRUE),
    max_cs_size = max(cs_size, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_tbl)

## =========================
## OPTIONAL: NULL CALIBRATION DIAGNOSTICS
## =========================

## Distribution of CS sizes under null
library(ggplot2)

ggplot(df_long %>% filter(!is.na(cs_size)),
       aes(x = cs_size, fill = method)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  scale_x_log10()+
  facet_grid(nullweight ~ thresh_lowcount) +
  #scale_x_continuous(trans = "log1p") +
  theme_bw() +
  labs(
    title = "FSuSiE credible set size under null",
    x = "CS size (log1p scale)",
    y = "Count"
  )

## Probability of detecting any CS under null
ggplot(df_long %>%
         distinct(method, nullweight, thresh_lowcount, has_cs),
       aes(x = factor(thresh_lowcount),
           y = has_cs,
           color = factor(nullweight),
           group = nullweight)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun = mean, geom = "point") +
  facet_wrap(~ method) +
  theme_bw() +
  labs(
    y = "P(any CS | null)",
    x = "thresh_lowcount",
    color = "nullweight"
  )
