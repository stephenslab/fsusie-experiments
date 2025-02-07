

rm(list = ls())
data = readRDS("Fig4_data.rds")
library("ComplexUpset")
library("tidyverse")
plot_df = data$a1
# Create the plot
a1 = ggplot(plot_df, aes(y = context)) +
  # Stacked bar for total CS
  geom_col(aes(x = n_of_cs), fill = "grey70", width = 0.6) +
  
  # Overlay the super_finemapped portion as a darker bar
  geom_col(aes(x = total_super_finemapped), fill = "steelblue", width = 0.6) +
  
  # Add text for total_super_finemapped inside the colored portion
  # We'll place it at half the super finemapped width.
  geom_text(aes(x = total_super_finemapped/2, 
                label = total_super_finemapped), 
            color = "white", size = 15, fontface = "bold") +
  
  # Add text for total number of CS to the right of the entire bar
  # We'll place it a bit to the right of the bar
  geom_text(aes(x = n_of_cs + (max(n_of_cs)*0.03), 
                label = n_of_cs), 
            color = "black", size = 15, fontface = "bold") +
  
  # Make it horizontal, adjust theme
  coord_flip() +
  
  # Add some spacing so text is not cut off
  scale_x_continuous(expand = expansion(mult = c(0,0.1))) +
  
  theme_minimal(base_size = 24) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 40),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) 
a1

custom_colors <- c(
  "cell-type specific" = "steelblue"
)
upset_df = data$a2
a2 = ComplexUpset::upset(
  upset_df%>%as_tibble,
  colnames(upset_df%>%as_tibble%>%select(-class,-cell_type)),
  keep_empty_groups = FALSE,
  base_annotations = list(
    `Intersection size` = intersection_size(
      mapping = aes(fill = as.character(upset_df$cell_type)),  # Reverse the order for bottom-to-top stacking
      bar_number_threshold = 1,
      width = 0.7,
      text = list(size = 10),
    ) + 
    ylab("") +ylim(c(0,13000))+ theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))+
    scale_fill_manual(na.value = "gray70" ,values = rev(custom_colors),name = "")
  ),
  set_sizes = upset_set_size(
    geom = geom_bar(aes(),fill = "gray70"), 
    filter_intersections = FALSE
  )+theme(axis.text.x = element_text(size = 10)),
  width_ratio = 0.15,
  themes = upset_default_themes(
    plot.margin = unit(c(0, 0, 0, 0), "mm"),
    axis.text = element_text(size = 30),
    axis.title.x = element_blank(),
    text = element_text(size = 40) ,        legend.position = "top"  # **Added this line**
 ),
  min_degree = 1
)+
  theme(legend.position = "top")  # Added this line       
a2
plot_df = data$b1
# Calculate Pearson correlation
cor_test <- cor.test(plot_df%>%pull(top_z.x), plot_df%>%pull(top_z.y))
cor_coeff <- round(cor_test$estimate, 2)
p_value <- signif(cor_test$p.value, 2)
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(cowplot)
library(RColorBrewer)
palette_overlap <- c("Overlap" = "steelblue", "No Overlap" = "grey70")

cor_coeff <- round(cor_test$estimate, 2)
p_value <- signif(cor_test$p.value, 2)

# Create the plot
b1 <- ggplot(plot_df, aes(x = top_z.x, y = top_z.y)) +
  geom_point(aes(),size = 10,color = "steelblue", size = 3, alpha = 0.8) +  # Adjusted size for readability
  geom_smooth(size = 3,color = "black",
    aes(),
    method = "lm",
    se = F,  # Set to TRUE to include confidence intervals
    size = 3  # Adjusted line size for clarity
  ) +
  scale_color_manual(
    values = palette_overlap,
    name = ""  # Updated legend title for clarity
  ) +
  scale_linetype_manual(
    values = c("eQTL" = "dashed", "Non eQTL" = "solid"),
    name = ""  # Added name for linetype legend
  ) +
  scale_shape_manual(
    values = c("eQTL" = 16, "Non eQTL" = 17),
    name = ""  # Added name for shape legend
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right"
  ) +
  labs(
    x = "Estimated Effect (car QTL)",
    y = "Estimated Effect (dmr QTL)"
  ) +
  annotate(
    "text",
    x = min(plot_df%>%pull(top_z.x))*1.2,
    y = max(plot_df%>%pull(top_z.y))*0.9,
    label = paste0("Pearson r = ", cor_coeff, "\nP-value = ", p_value),
    hjust = 0,
    size = 10,
    color = "black"
  )

# Display the plot
print(b1)

plot_df = data$b2
# Calculate Pearson correlation
cor_test <- cor.test(plot_df%>%pull(top_z.x), plot_df%>%pull(top_z.y))
cor_coeff <- round(cor_test$estimate, 2)
p_value <- signif(cor_test$p.value, 2)
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(cowplot)
library(RColorBrewer)
palette_overlap <- c("Overlap" = "steelblue", "No Overlap" = "grey70")

cor_coeff <- round(cor_test$estimate, 2)
p_value <- signif(cor_test$p.value, 2)

# Create the plot
b2 <- ggplot(plot_df, aes(x = top_z.x, y = top_z.y)) +
  geom_point(aes(),size = 10,color = "grey", size = 3, alpha = 0.8) +  # Adjusted size for readability
  geom_smooth(size = 3,color = "black",
    aes(),
    method = "lm",
    se = F,  # Set to TRUE to include confidence intervals
    size = 3  # Adjusted line size for clarity
  ) +
  scale_color_manual(
    values = palette_overlap,
    name = ""  # Updated legend title for clarity
  ) +
  scale_linetype_manual(
    values = c("eQTL" = "dashed", "Non eQTL" = "solid"),
    name = ""  # Added name for linetype legend
  ) +
  scale_shape_manual(
    values = c("eQTL" = 16, "Non eQTL" = 17),
    name = ""  # Added name for shape legend
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right"
  ) +
  labs(
    x = "Estimated Effect (car QTL)",
    y = "Estimated Effect (dmr QTL)"
  ) +
  annotate(
    "text",
    x = min(plot_df%>%pull(top_z.x))*1.2,
    y = max(plot_df%>%pull(top_z.y))*0.9,
    label = paste0("Pearson r = ", cor_coeff, "\nP-value = ", p_value),
    hjust = 0,
    size = 10,
    color = "black"
  )

# Display the plot
print(b2)

plot_df = data$b3
# Calculate Pearson correlation
cor_test <- cor.test(plot_df%>%pull(epi_top_z), plot_df%>%pull(betahat))
cor_coeff <- round(cor_test$estimate, 2)
p_value <- signif(cor_test$p.value, 2)
b3 <- plot_df%>%
  ggplot(aes(y = epi_top_z, x = betahat)) +
  geom_point( size = 10,color = "steelblue", alpha = 0.5) +
  geom_smooth(aes( ), method = "lm", se = TRUE, size = 3, color = "black") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right"
  ) +
  labs(
    y = "Estimated Effect (dmr QTL)",
    x = "Estimated Effect (eQTL)"
  ) +
  annotate(
    "text",
    x = min(plot_df$epi_top_z, na.rm = TRUE) + 0.2,
    y = max(plot_df %>% filter(context.x =="ROSMAP_DLPFC_mQTL" )%>%pull(epi_top_z), na.rm = TRUE) - 0.05,
    label = paste0("Pearson r = ", cor_coeff, "\nP-value = ", p_value),
    hjust = 0,
    size = 10
  )+ylim(-0.5,0.5)
# Calculate Pearson correlation
b3

plot_df = data$b4
# Calculate Pearson correlation
cor_test <- cor.test(plot_df%>%pull(epi_top_z), plot_df%>%pull(betahat))
cor_coeff <- round(cor_test$estimate, 2)
p_value <- signif(cor_test$p.value, 2)
b4 <- plot_df%>%
  ggplot(aes(y = epi_top_z, x = betahat)) +
  geom_point( size = 10,color = "maroon", alpha = 0.5) +
  geom_smooth(aes( ), method = "lm", se = TRUE, size = 3, color = "black") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right"
  ) +
  labs(
    y = "Top Estimated Effect (car QTL)",
    x = "Median Estimated Effect (eQTL)"
  ) +
  annotate(
    "text",
    x = min(plot_df$epi_top_z, na.rm = TRUE) + 0.2,
    y = max(plot_df$betahat, na.rm = TRUE) - 0.05,
    label = paste0("Pearson r = ", cor_coeff, "\nP-value = ", p_value),
    hjust = 0,
    size = 10
  )
b4

# Interface: Assign datasets dynamically

plot_df <- data$d[[1]]
haQTL_df <- data$d[[2]]
mQTL_df <- data$d[[3]]
gene_info <- data$d[[4]]
sumstat <- data$d[[5]]
pip_df <- data$d[[6]]

# Parameters
view_win <- c(4759843, 5000000)
text_size <- 20
gene <- c("ENSG00000029725", "ENSG00000161929", "ENSG00000108556")

custom_labeller <- function(x) {
  x %>% 
    gsub("DeJager_", "", ., fixed = TRUE) %>%  
    gsub("([_:,|-])", "\n", .)             
}

# Plot 1
p1 <- ggplot() +
  geom_point(data = plot_df %>% filter(study == "DLPFC_DeJager_eQTL"), 
             aes(x = pos, y = `z`), size = 8, alpha = 0.1) +
  facet_grid(study + region ~ ., labeller = labeller(.rows = custom_labeller), scale = "free_y") +
  geom_point(data = plot_df %>% filter(study == "DLPFC_DeJager_eQTL", CS1),
             aes(x = pos, y = `z`), color = "#2E8B57", size = 8, alpha = 1) +
  geom_point(data = plot_df %>% filter(CS4, study == "DLPFC_DeJager_eQTL"),
             aes(x = pos, y = `z`), color = "steelblue", size = 8, alpha = 1) +
  geom_point(data = plot_df %>% filter(CS1, study == "DLPFC_DeJager_eQTL"),
             aes(x = pos, y = `z`), color = "maroon", size = 8, alpha = 1) +
  geom_point(data = plot_df %>% filter(variant_alternate_id %in% (plot_df %>% filter(CS2, study == "Mic_mega_eQTL") %>% pull(variant_alternate_id)), study == "AD_Bellenguez_2022"),
             aes(x = pos, y = `z`), color = "#2E8B57", size = 8, alpha = 1, shape = "triangle") +
  xlim(view_win) +
  theme_bw() +
  theme(text = element_text(size = text_size),
        strip.text.y = element_text(size = text_size, angle = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = text_size),
        axis.title = element_text(size = text_size)) +
  xlab("") +
  ylab("z")

# Plot 2
p2 <- ggplot() + theme_bw() +
  xlim(view_win) +
  ylab("Estimated effect") +
  geom_line(data = haQTL_df %>% mutate(study = "haQTL effect") %>% filter(CS == 2),
            aes_string(y = "fun_plot", x = "x", col = "CS"), size = 4, col = "steelblue") +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(aes(x = view_win + 1000, xend = (haQTL_df %>% mutate(study = "haQTL effect") %>% filter(CS == 2))$start[[1]], y = 0, yend = 0, col = "CS"), size = 4, col = "steelblue") +
  geom_segment(aes(xend = view_win - 1000, x = (haQTL_df %>% mutate(study = "haQTL effect") %>% filter(CS == 2))$end[[1]], y = 0, yend = 0, col = "CS"), size = 4, col = "steelblue") +
  geom_line(data = mQTL_df %>% mutate(study = "haQTL effect") %>% filter(CS == 7),
            aes_string(y = "fun_plot", x = "x", col = "CS"), size = 4, col = "maroon") +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(aes(x = view_win + 1000, xend = (mQTL_df %>% mutate(study = "haQTL effect") %>% filter(CS == 7))$start[[1]], y = 0, yend = 0, col = "CS"), size = 4, col = "maroon") +
  geom_segment(aes(xend = view_win - 1000, x = (mQTL_df %>% mutate(study = "haQTL effect") %>% filter(CS == 7))$end[[1]], y = 0, yend = 0, col = "CS"), size = 4, col = "maroon") +
  geom_hline(aes(yintercept = 0)) +
  theme(text = element_text(size = text_size), 
        strip.text.y = element_text(size = text_size, angle = 0.5),
        axis.text.y = element_text(size = text_size),
        axis.title.x = element_text(size = text_size)) +
  xlab("Phenotype Position") +
  ylab("Estimated\neffect") +
  geom_segment(arrow = arrow(length = unit(0.5, "cm")), aes(x = start, xend = end, y = -0.3, yend = -0.3), size = 1,
               data = gene_info %>% filter(gene_id %in% c(gene)) %>% mutate(study = "gene_plot")) +
  geom_text(aes(x = (start + end) / 2, y = -0.25, label = gene_name), size = 10, 
            data = gene_info %>% filter(gene_id %in% gene) %>% mutate(study = "gene_plot")) + 
  geom_point(data = sumstat %>% filter(ha), aes(x = pos - start_distance, y = beta), color = "steelblue", size = 2) +
  geom_point(data = sumstat %>% filter(!ha), aes(x = pos - start_distance, y = beta), color = "maroon", size = 2) +
  ylim(c(-0.31, 0.7))

# Plot 3
p3 <- ggplot() +
  facet_grid(study ~ ., scale = "fixed", labeller = labeller(study = c(ROSMAP_DLPFC_haQTL = "car-QTL", ROSMAP_DLPFC_mQTL = "dmr-QTL"))) +
  xlim(view_win) +
  theme_bw() +
  ylim(c(0, 0.4)) +
  theme(text = element_text(size = text_size),
        strip.text.y = element_text(size = text_size, angle = 0.5),
        axis.text.y = element_text(size = text_size),
        axis.title.x = element_text(size = text_size)) +
  scale_y_continuous(breaks = pretty((pip_df %>% filter(study %in% c("ROSMAP_DLPFC_haQTL", "ROSMAP_DLPFC_mQTL", ""), cs_coverage_0.95_min_corr == 2))$pip, n = 3)) +
  xlab("Genotype Position") +
  ylab("PIP") +
  geom_point(data = pip_df %>% filter(study %in% c("ROSMAP_DLPFC_haQTL"), cs_coverage_0.95_min_corr == 2),
             aes(x = pos, y = pip, color = as.character(cs_coverage_0.95)), alpha = 1, size = 8, color = "steelblue") +
  geom_point(data = pip_df %>% filter(study %in% c("ROSMAP_DLPFC_mQTL", ""), cs_coverage_0.95 == 7),
             aes(x = pos, y = pip, color = as.character(cs_coverage_0.95)), alpha = 1, size = 14, color = "maroon")

d <- cowplot::plot_grid(plotlist = list(p1, p3, p2), ncol = 1, align = "v", axis = "tlbr", label_size = 45, label_fontface = "bold", rel_heights = c(4, 3, 4))
d



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
# Plot 1
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
# Plot 2
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
# Plot 3
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

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)      
library(magrittr)   

# Interface: Assign datasets dynamically

plot_df <- data$f[[1]]
haQTL_df <- data$f[[2]]
MSBB_df <- data$f[[3]]
gene_info <- data$f[[4]]
sumstat <- data$f[[5]]
pip_df <- data$f[[6]]
QTL_data <- data$f[[7]] 
# Custom Labeller
custom_labeller <- function(x) {
  x %>% 
    gsub("mega_", "", ., fixed = TRUE) %>%  
    gsub("([_:,|-])", " \n ", .)             
}

# Parameters
text_size <- 20
view_win <- c(5.12e7, 5.16e7)
gene <- c("ENSG00000139629", "ENSG00000050438")

# -----------------
# Plot 1
# -----------------
p1 <- ggplot() +
  geom_point(
    data = plot_df %>% filter(),
    aes(x = pos, y = `z`),
    size = 8, alpha = 0.1
  ) +
  facet_grid(
    study + region ~ .,
    labeller = labeller(.rows =custom_labeller),
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
  ylab("Z")

# -----------------
# Plot 2
# -----------------
p2 <- ggplot() +
  theme_bw() +
  theme(
    text = element_text(size = text_size),
    strip.text.y = element_text(size = text_size, angle = 0.5),
    axis.text.x = element_text(size = text_size),    
    axis.text.y = element_text(size = text_size),
    axis.title.x = element_text(size = text_size)
  ) +
  xlim(view_win) +
  geom_line(
    data = haQTL_df %>% mutate(study = "haQTL effect") %>% filter(CS == 5),
    aes_string(y = "fun_plot", x = "x", col = "CS"),
    size = 4, col = "steelblue"
  ) +
  geom_line(
    data = MSBB_df %>% mutate(study = "dmr QTL effect") %>% filter(CS == 14),
    aes_string(y = "fun_plot", x = "x", col = "CS"),
    size = 4, col = "maroon"
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
  ylab("Estimated\neffect")  +
  geom_point(data =  sumstat %>% filter(ha), aes(x = pos - start_distance, y = beta), color = "steelblue", size = 2) +
  geom_point(data =  sumstat %>% filter(!ha), aes(x = pos - start_distance, y = beta), color = "maroon", size = 2)

# -----------------
# Plot 3
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
    axis.text.x = element_text(size = text_size),
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
  geom_point(
    data = QTL_data %>%
      filter(variant_id == "chr12:51362485:T:C", study == "MSBB_mQTL") %>%
      mutate(study = "dmr-QTL") %>%
      separate(col = variant_id, into = c("chrom", "pos"), remove = FALSE) %>%
      mutate(pos = as.numeric(pos)),
    aes(x = pos, y = pip, color = as.character(cs_coverage_0.95)),
    alpha = 1, size = 8, color = "maroon"
  ) +
  xlab("Genotype Position") +
  ylab("PIP") 

# Combine Plots
f <- cowplot::plot_grid(plotlist = list(p1, p3, p2),
  ncol = 1,
  align = "v",
  axis = "tlbr",
  label_fontface = "bold",
  rel_heights = c(5, 2, 4)
) 

f

merged_A = plot_grid(a1, a2+theme(text = element_text(15)), ncol = 2, align = 'h', labels = NULL,rel_widths = c(0.3, 0.7) )
merged_B = cowplot::plot_grid(b1,b2,b3,b4, nrow = 2,ncol = 2) +
  theme_cowplot(font_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5),  # Increased title size for emphasis
    axis.text = element_text( size = 24),
        axis.title = element_text( size = 24),
    legend.title = element_text( size = 24),
    legend.text = element_text(size = 24),
    legend.position = "top"
  )
ab = plot_grid( merged_A , merged_B+theme(plot.margin = unit(c(0, 0,0,1.5), "cm"))
 , labels = c("A", "B") ,label_size = 80,rel_widths = c(1, 1))
cd  = plot_grid( d+theme(plot.margin = unit(c(0, 0,0,1.5), "cm")) , d+theme(plot.margin = unit(c(0, 0,0,1.5), "cm"))
 , labels = c("C", "D") ,label_size = 80,rel_widths = c(1, 1))
ef  = plot_grid( e+theme(plot.margin = unit(c(0, 0,0,1.5), "cm")) , f+theme(plot.margin = unit(c(0, 0,0,1.5), "cm"))
 , labels = c("E", "F") ,label_size = 80,rel_widths = c(1, 1))

options(repr.plot.width = 50, repr.plot.height = 45)  

plot_grid( ab , cd,ef
 , rel_widths = c(1, 1), nrow = 3)

