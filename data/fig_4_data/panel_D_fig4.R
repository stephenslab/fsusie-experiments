
# Interface: Assign datasets dynamically
plot_df <- data$e[[1]]
haQTL_df <- data$e[[2]]
gene_info <- data$e[[3]]
sumstat <- data$e[[4]]
pip_df <- data$e[[5]]

# Custom Labeller
custom_labeller <- function(x) {
  x %>% 
    gsub("DeJager_", "", ., fixed = TRUE) %>%  
    gsub("([_:,|-])", " \n ", .)             
}

# Parameters
text_size <- 20
view_win <- c(207317782, 207895513)
gene <- c("chr1:205117782-208795513", "1_205972031_208461272", "ENSG00000203710","ENSG00000117322")

# -----------------
# Plot 1  ------
# -----------------
p1 <- ggplot() +
  geom_point(
    data = plot_df %>% filter(),
    aes(x = pos, y = `z`),
    size = 8, alpha = 0.1
  ) +
  facet_grid(
    study + region ~ .,
    labeller = labeller(.rows = custom_labeller),
    scale = "free_y"
  ) +
  geom_point(
    data = plot_df %>% filter(CS1),
    aes(x = pos, y = `z`),
    color = "steelblue",
    size = 8, alpha = 1
  ) +
  xlim(view_win) +
  theme_bw() +
  theme(
    text = element_text(size = text_size),
    strip.text.y = element_text(size = text_size, angle = 0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = text_size),
    axis.title = element_text(size = text_size)
  ) +
  xlab("") +
  ylab("z")

# -----------------
# Plot 2------
# -----------------
p2 <- ggplot() +
  theme_bw() +
  theme(
    text = element_text(size = text_size),
    strip.text.y = element_text(size = text_size, angle = 0.5),
    axis.text.x = element_text(size = text_size),
    axis.title.x = element_text(size = text_size)
  ) +
  xlim(view_win) +
  ylab("Estimated effect") +
  geom_line(
    data = haQTL_df %>% mutate(study = "haQTL effect") %>% filter(CS == 5),
    aes_string(y = "fun_plot", x = "x", col = "CS"),
    size = 4, col = "steelblue"
  ) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(
    aes(x = view_win[2] + 1000, 
        xend = (haQTL_df %>% mutate(study = "haQTL effect") %>% filter(CS == 5))$start[[1]],
        y = 0, yend = 0, col = "CS"),
    size = 4, col = "steelblue"
  ) +
  geom_segment(
    aes(xend = view_win[1] - 1000,
        x = (haQTL_df %>% mutate(study = "haQTL effect") %>% filter(CS == 5))$end[[1]],
        y = 0, yend = 0, col = "CS"),
    size = 4, col = "steelblue"
  ) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(
    arrow = arrow(length = unit(0.5, "cm")), 
    aes(x = start, xend = end, y = -0.12 - strand / 10 - 0.02, 
        yend = -0.12 - strand / 10 - 0.02),
    size = 1,
    data = gene_info %>%
      filter(gene_id %in% gene) %>%
      mutate(study = "gene_plot")
  ) +
  geom_text(
    aes(x = (start + end) / 2, 
        y = -0.12 - strand / 10,
        label = gene_name),
    size = 10, 
    data = gene_info %>%
      filter(gene_id %in% gene) %>%
      mutate(study = "gene_plot")
  ) +
  xlab("Phenotype Position") +
  ylab("Estimated\neffect") +
  geom_point(data = sumstat, aes(x = pos - start_distance, y = beta), color = "steelblue", size = 2)

# -----------------
# Plot 3 -----
# -----------------
p3 <- ggplot() +
  facet_grid(
    study ~ .,
    scale = "free_y",
    labeller = labeller(study = c(ROSMAP_DLPFC_haQTL = "car-QTL"))
  ) +
  xlim(view_win) +
  theme_bw() +
  theme(
    text = element_text(size = text_size),
    strip.text.y = element_text(size = text_size, angle = 0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = text_size),
    axis.title = element_text(size = text_size)
  ) +
  scale_y_continuous(
    breaks = pretty(c(0, 1), n = 3),
    limits = c(0, 1)
  ) +
  xlab("") +
  ylab("PIP") +
  geom_point(
    data = pip_df %>%
      filter(study %in% c("ROSMAP_DLPFC_haQTL"), cs_coverage_0.95 == 5),
    aes(x = pos, y = pip, color = as.character(cs_coverage_0.95)),
    alpha = 1, size = 8, color = "steelblue"
  ) +
  xlab("Genotype Position") +
  ylab("PIP") 

# Combine Plots
e <- cowplot::plot_grid(plotlist = list(p1, p3, p2),
                        ncol = 1,
                        align = "v",
                        axis = "tlbr",
                        label_size = 45,  
                        label_fontface = "bold",
                        rel_heights = c(4, 2, 4)
) 

e 
