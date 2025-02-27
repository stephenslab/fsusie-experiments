rm(list = ls())

library(ComplexUpset)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(tidyr)
library(magrittr)

# Define common text size settings for consistencyhttp://127.0.0.1:15077/graphics/4ebeea29-43a3-49b0-bdce-8e9f155ba977.png
axis_text_size <- 10      # for axis texts and theme elements
annotation_size <- 5      # for numeric labels on bars

# Set working directory and load data
path <- getwd()
data <- readRDS(file.path(path, "data/fig_4_data/Fig4_data.rds"))

### Plot a1: Horizontal Stacked Bar Plot ###
plot_df <- data$a1

plot_df$context[1] = "car QTL"
plot_df$context[2] = "dmr QTL"

a1 <-ggplot(plot_df, aes(y = context)) +
  # Draw the full-length bar (total number of CS)
  geom_col(aes(x = n_of_cs, fill = " "), width = 0.6) +
  # Overlay the portion that is super-finemapped
  geom_col(aes(x = total_super_finemapped, fill = "1 SNP CS"), width = 0.6) +
  # Label the super-finemapped portion (centered inside the green bar)
  geom_text(aes(x = total_super_finemapped/2, label = total_super_finemapped),
            color = "white", size = annotation_size, fontface = "bold") +
  # Label the total CS count (positioned slightly to the right of the bar)
  geom_text(aes(x = n_of_cs + (max(n_of_cs)*0.03), label = n_of_cs),
            color = "black", size = annotation_size, fontface = "bold") +
  coord_flip() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c("Total CS" = "grey70", "1 SNP CS" = "green4")) +
  theme_minimal(base_size = axis_text_size) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = axis_text_size),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 10, 10, 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.text.x.bottom =  element_text(size = 14),
        legend.position = "bottom")  # Adjust as needed


### Plot a2: Upset Plot ###
custom_colors <- c("cell-type specific" = "royalblue")
upset_df <- data$a2

a2 <- ComplexUpset::upset(
  upset_df %>% as_tibble(),
  colnames(upset_df %>% as_tibble() %>% select(-class, -cell_type)),
  keep_empty_groups = FALSE,
  set_sizes = FALSE,
  base_annotations = list(
    `Intersection size` = intersection_size(
      mapping = aes(fill = as.character(upset_df$cell_type)),
      bar_number_threshold = 1,
      width = 0.7,
      text = list(size = annotation_size)
    ) +
      ylab("") +
      ylim(c(0, 13000)) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "mm")) +
      scale_fill_manual(na.value = "gray70",
                        values = rev(custom_colors),
                        name = "")
  ),
  width_ratio = 0.15,  # Adjust if needed
  themes = upset_default_themes(
    plot.margin = unit(c(0, 0, 0, 0), "mm"),
    axis.text = element_text(size = axis_text_size),
    axis.title.x = element_blank(),
    text = element_text(size = axis_text_size),
    legend.position = "bottom",
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 14)
  ),
  min_degree = 1
) +
  theme(legend.position = "bottom")


# Optionally, combine the two plots using cowplot for a unified display
combined_plot <- cowplot::plot_grid(a1, a2, ncol = 2, 
                                    rel_widths = c(0.2,1),
                                    rel_heights = c(0.2, 1)
                                    )
print(combined_plot)

save_path=  paste0(getwd(),
                   "/plot/"
)
ggsave(combined_plot, file=paste0(save_path,"upset.pdf"),
       width = 50,
       height =25,
       units = "cm"
)
#### Correlation plot ----

plot_df = data$b1
# Calculate Pearson correlation
cor_test <- cor.test(plot_df%>%pull(top_z.x), plot_df%>%pull(top_z.y))
cor_coeff <- round(cor_test$estimate, 2)
p_value <- signif(cor_test$p.value, 2)
# Load necessary libraries



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
    x = -1.5,
    y = 0.45,
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
P0= cowplot::plot_grid(b1,b2,b3,b4, nrow = 2,ncol = 2) 
save_path=  paste0(getwd(),
                   "/plot/"
)
ggsave(P0, file=paste0(save_path,"sup_correlation_effect.pdf"),
       width = 60,
       height =45,
       units = "cm"
)
