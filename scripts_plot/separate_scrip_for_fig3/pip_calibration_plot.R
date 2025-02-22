# Load necessary library
library(ggplot2)
library(dplyr)

save_path=  paste0(getwd(),
                   "/plot/fig3_separate_panel/"
)
path <- getwd()

## Gaussian ----
load(paste( path,"/simulation/Simulation_results/comparison_susie_fusie_128_sd1.RData", 
            sep=""))


true_lab <- do.call( c,
                     lapply(1: length(res),
                            
                            function( i) {
                              
                              a <-  rep( 0,   length(res[[i]]$susiF_pip))
                              a[res[[i]]$true_pos] <- 1
                              return(a)
                            }
                            
                     )
)  





score_fsusie <-  do.call( c, lapply( 1: length(res),
                                     function( i) res[[i]]$susiF_pip))
score_sp_fsusie <-  do.call( c, lapply( 1: length(res),
                                        function( i) res[[i]]$susiF_sp_pip))

score_susie <-  do.call( c, lapply( 1: length(res),
                                    function( i) res[[i]]$susie_rss_pip))



pip_scores <- score_fsusie 
# 
causal_flags <-true_lab

num_bins <- 11
bin_edges <- seq(0, 1, length.out = num_bins + 1)
bin_indices <- cut(pip_scores, breaks = bin_edges, include.lowest = TRUE, labels = FALSE)

# Calculate bin centers
bin_centers <- (bin_edges[-length(bin_edges)] + bin_edges[-1]) / 2

# Create a data frame including the bin index
pip_df <- data.frame(pip_scores, causal_flags, bin_indices)

# Compute summary statistics per bin and add bin centers
summary_df <- pip_df %>%
  group_by(bin_indices) %>%
  summarise(
    mean_pip = mean(pip_scores),
    eip = mean(causal_flags),
    count = n()
  ) %>%
  mutate(
    se = sqrt(eip * (1 - eip) / count),
    ci_lower = pmax(eip - 1.96 * se, 0),
    ci_upper = pmin(eip + 1.96 * se, 1),
    # Map the bin index to the corresponding bin center
    bin_center = bin_centers[as.numeric(bin_indices)]
  )

# Optionally remove bins with too few observations
#summary_df <- summary_df[summary_df$count >= 60, ]

# Plot using bin centers as the x-axis
summary_df= summary_df[-which(summary_df$count<60),]
# Plot PIP Calibration with Confidence Intervals
P11 =ggplot(summary_df, aes(x = mean_pip, y = eip)) +
  geom_point(color = "#D41159", size = 3) +  # Scatter plot of PIP vs. EIP
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.02, color = "#D41159") +  # CI bars
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Ideal calibration line
  labs( 
    x = "Mean PIP in Bin", y = "Empirical Inclusion Probability  ") +
  # ggtitle("fsusie normal mixture")+
  theme_minimal()




pip_scores <- score_sp_fsusie
# 
causal_flags <-true_lab


bin_edges <- seq(0, 1, length.out = num_bins + 1)
bin_indices <- cut(pip_scores, breaks = bin_edges, include.lowest = TRUE, labels = FALSE)

# Calculate bin centers
bin_centers <- (bin_edges[-length(bin_edges)] + bin_edges[-1]) / 2

# Create a data frame including the bin index
pip_df <- data.frame(pip_scores, causal_flags, bin_indices)

# Compute summary statistics per bin and add bin centers
summary_df <- pip_df %>%
  group_by(bin_indices) %>%
  summarise(
    mean_pip = mean(pip_scores),
    eip = mean(causal_flags),
    count = n()
  ) %>%
  mutate(
    se = sqrt(eip * (1 - eip) / count),
    ci_lower = pmax(eip - 1.96 * se, 0),
    ci_upper = pmin(eip + 1.96 * se, 1),
    # Map the bin index to the corresponding bin center
    bin_center = bin_centers[as.numeric(bin_indices)]
  )

# Optionally remove bins with too few observations
#summary_df <- summary_df[summary_df$count >= 60, ]

# Plot using bin centers as the x-axis
summary_df= summary_df[-which(summary_df$count<60),]
# Plot PIP Calibration with Confidence Intervals
 


 
# Plot PIP Calibration with Confidence Intervals
P21 =ggplot(summary_df, aes(x = bin_center, y = eip)) +
  geom_point(color = "#1A85FF", size = 3) +  # Scatter plot of PIP vs. EIP
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.02, color = "#1A85FF") +  # CI bars
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Ideal calibration line
  labs( 
    x = "Mean PIP in Bin", y = "Empirical Inclusion Probability  ") +
  # ggtitle("fsusie  mixture per scale")+
  theme_minimal()






pip_scores <-score_susie
# 
causal_flags <-true_lab

# Define bins

bin_edges <- seq(0, 1, length.out = num_bins + 1)
bin_indices <- cut(pip_scores, breaks = bin_edges, include.lowest = TRUE, labels = FALSE)

# Calculate bin centers
bin_centers <- (bin_edges[-length(bin_edges)] + bin_edges[-1]) / 2

# Create a data frame including the bin index
pip_df <- data.frame(pip_scores, causal_flags, bin_indices)

# Compute summary statistics per bin and add bin centers
summary_df <- pip_df %>%
  group_by(bin_indices) %>%
  summarise(
    mean_pip = mean(pip_scores),
    eip = mean(causal_flags),
    count = n()
  ) %>%
  mutate(
    se = sqrt(eip * (1 - eip) / count),
    ci_lower = pmax(eip - 1.96 * se, 0),
    ci_upper = pmin(eip + 1.96 * se, 1),
    # Map the bin index to the corresponding bin center
    bin_center = bin_centers[as.numeric(bin_indices)]
  )

# Optionally remove bins with too few observations
#summary_df <- summary_df[summary_df$count >= 60, ]

# Plot using bin centers as the x-axis
 

summary_df= summary_df[-which(summary_df$count<60),]
# Plot PIP Calibration with Confidence Intervals
P31 =ggplot(summary_df, aes(x = bin_center, y = eip)) +
  geom_point(color = "#40B0A6", size = 3) +  # Scatter plot of PIP vs. EIP
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.02, color = "#40B0A6") +  # CI bars
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Ideal calibration line
  labs( 
    x = "Mean PIP in Bin", y = "Empirical Inclusion Probability  ") +
  # ggtitle("fsusie normal mixture")+
  theme_minimal()







##### Block ----
path <- getwd()
load(paste( path,"/simulation/Simulation_results/comparison_susie_fusie_block_sd1.RData", 
            sep=""))


true_lab <- do.call( c,
                     lapply(1: length(res),
                            
                            function( i) {
                              
                              a <-  rep( 0,   length(res[[i]]$susiF_pip))
                              a[res[[i]]$true_pos] <- 1
                              return(a)
                            }
                            
                     )
)  





score_fsusie <-  do.call( c, lapply( 1: length(res),
                                     function( i) res[[i]]$susiF_pip))
score_sp_fsusie <-  do.call( c, lapply( 1: length(res),
                                        function( i) res[[i]]$susiF_sp_pip))

score_susie <-  do.call( c, lapply( 1: length(res),
                                    function( i) res[[i]]$susie_rss_pip))



pip_scores <- score_fsusie 
# 
causal_flags <-true_lab

# Define bins

bin_edges <- seq(0, 1, length.out = num_bins + 1)
bin_indices <- cut(pip_scores, breaks = bin_edges, include.lowest = TRUE, labels = FALSE)

# Calculate bin centers
bin_centers <- (bin_edges[-length(bin_edges)] + bin_edges[-1]) / 2

# Create a data frame including the bin index
pip_df <- data.frame(pip_scores, causal_flags, bin_indices)

# Compute summary statistics per bin and add bin centers
summary_df <- pip_df %>%
  group_by(bin_indices) %>%
  summarise(
    mean_pip = mean(pip_scores),
    eip = mean(causal_flags),
    count = n()
  ) %>%
  mutate(
    se = sqrt(eip * (1 - eip) / count),
    ci_lower = pmax(eip - 1.96 * se, 0),
    ci_upper = pmin(eip + 1.96 * se, 1),
    # Map the bin index to the corresponding bin center
    bin_center = bin_centers[as.numeric(bin_indices)]
  )

# Optionally remove bins with too few observations
#summary_df <- summary_df[summary_df$count >= 60, ]

 


summary_df= summary_df[-which(summary_df$count<60),]
# Plot PIP Calibration with Confidence Intervals
P12 =ggplot(summary_df, aes(x =bin_center, y = eip)) +
  geom_point(color = "#D41159", size = 3) +  # Scatter plot of PIP vs. EIP
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.02, color = "#D41159") +  # CI bars
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Ideal calibration line
  labs( 
    x = "Mean PIP in Bin", y = "Empirical Inclusion Probability  ") +
  # ggtitle("fsusie normal mixture")+
  theme_minimal()




pip_scores <- score_sp_fsusie
# 
causal_flags <-true_lab

# Define bins
# Define bins

bin_edges <- seq(0, 1, length.out = num_bins + 1)
bin_indices <- cut(pip_scores, breaks = bin_edges, include.lowest = TRUE, labels = FALSE)

# Calculate bin centers
bin_centers <- (bin_edges[-length(bin_edges)] + bin_edges[-1]) / 2

# Create a data frame including the bin index
pip_df <- data.frame(pip_scores, causal_flags, bin_indices)

# Compute summary statistics per bin and add bin centers
summary_df <- pip_df %>%
  group_by(bin_indices) %>%
  summarise(
    mean_pip = mean(pip_scores),
    eip = mean(causal_flags),
    count = n()
  ) %>%
  mutate(
    se = sqrt(eip * (1 - eip) / count),
    ci_lower = pmax(eip - 1.96 * se, 0),
    ci_upper = pmin(eip + 1.96 * se, 1),
    # Map the bin index to the corresponding bin center
    bin_center = bin_centers[as.numeric(bin_indices)]
  )



summary_df= summary_df[-which(summary_df$count<60),]
# Plot PIP Calibration with Confidence Intervals
P22 =ggplot(summary_df, aes(x = bin_center, y = eip)) +
  geom_point(color = "#1A85FF", size = 3) +  # Scatter plot of PIP vs. EIP
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.02, color = "#1A85FF") +  # CI bars
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Ideal calibration line
  labs( 
    x = "Mean PIP in Bin", y = "Empirical Inclusion Probability  ") +
  # ggtitle("fsusie  mixture per scale")+
  theme_minimal()






pip_scores <-score_susie
# 
causal_flags <-true_lab

# Define bins

bin_edges <- seq(0, 1, length.out = num_bins + 1)
bin_indices <- cut(pip_scores, breaks = bin_edges, include.lowest = TRUE, labels = FALSE)

# Calculate bin centers
bin_centers <- (bin_edges[-length(bin_edges)] + bin_edges[-1]) / 2

# Create a data frame including the bin index
pip_df <- data.frame(pip_scores, causal_flags, bin_indices)

# Compute summary statistics per bin and add bin centers
summary_df <- pip_df %>%
  group_by(bin_indices) %>%
  summarise(
    mean_pip = mean(pip_scores),
    eip = mean(causal_flags),
    count = n()
  ) %>%
  mutate(
    se = sqrt(eip * (1 - eip) / count),
    ci_lower = pmax(eip - 1.96 * se, 0),
    ci_upper = pmin(eip + 1.96 * se, 1),
    # Map the bin index to the corresponding bin center
    bin_center = bin_centers[as.numeric(bin_indices)]
  )


summary_df= summary_df[-which(summary_df$count<60),]
# Plot PIP Calibration with Confidence Intervals
P32 =ggplot(summary_df, aes(x = bin_center, y = eip)) +
  geom_point(color = "#40B0A6", size = 3) +  # Scatter plot of PIP vs. EIP
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.02, color = "#40B0A6") +  # CI bars
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Ideal calibration line
  labs( 
    x = "Mean PIP in Bin", y = "Empirical Inclusion Probability  ") +
  # ggtitle("fsusie normal mixture")+
  theme_minimal()






#### Decay ----


path <- getwd()
load(paste( path,"/simulation/Simulation_results/comparison_susie_fusie_distdecay_sd1.RData", 
            sep=""))


true_lab <- do.call( c,
                     lapply(1: length(res),
                            
                            function( i) {
                              
                              a <-  rep( 0,   length(res[[i]]$susiF_pip))
                              a[res[[i]]$true_pos] <- 1
                              return(a)
                            }
                            
                     )
)  





score_fsusie <-  do.call( c, lapply( 1: length(res),
                                     function( i) res[[i]]$susiF_pip))
score_sp_fsusie <-  do.call( c, lapply( 1: length(res),
                                        function( i) res[[i]]$susiF_sp_pip))

score_susie <-  do.call( c, lapply( 1: length(res),
                                    function( i) res[[i]]$susie_rss_pip))



pip_scores <- score_fsusie 
# 
causal_flags <-true_lab

# Define bins

bin_edges <- seq(0, 1, length.out = num_bins + 1)
bin_indices <- cut(pip_scores, breaks = bin_edges, include.lowest = TRUE, labels = FALSE)

# Calculate bin centers
bin_centers <- (bin_edges[-length(bin_edges)] + bin_edges[-1]) / 2

# Create a data frame including the bin index
pip_df <- data.frame(pip_scores, causal_flags, bin_indices)

# Compute summary statistics per bin and add bin centers
summary_df <- pip_df %>%
  group_by(bin_indices) %>%
  summarise(
    mean_pip = mean(pip_scores),
    eip = mean(causal_flags),
    count = n()
  ) %>%
  mutate(
    se = sqrt(eip * (1 - eip) / count),
    ci_lower = pmax(eip - 1.96 * se, 0),
    ci_upper = pmin(eip + 1.96 * se, 1),
    # Map the bin index to the corresponding bin center
    bin_center = bin_centers[as.numeric(bin_indices)]
  )

summary_df= summary_df[-which(summary_df$count<60),]
# Plot PIP Calibration with Confidence Intervals
P13 =ggplot(summary_df, aes(x = bin_center, y = eip)) +
  geom_point(color = "#D41159", size = 3) +  # Scatter plot of PIP vs. EIP
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.02, color = "#D41159") +  # CI bars
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Ideal calibration line
  labs( 
    x = "Mean PIP in Bin", y = "Empirical Inclusion Probability  ") +
  # ggtitle("fsusie normal mixture")+
  theme_minimal()




pip_scores <- score_sp_fsusie
# 
causal_flags <-true_lab

# Define bins
# Define bins

bin_edges <- seq(0, 1, length.out = num_bins + 1)
bin_indices <- cut(pip_scores, breaks = bin_edges, include.lowest = TRUE, labels = FALSE)

# Calculate bin centers
bin_centers <- (bin_edges[-length(bin_edges)] + bin_edges[-1]) / 2

# Create a data frame including the bin index
pip_df <- data.frame(pip_scores, causal_flags, bin_indices)

# Compute summary statistics per bin and add bin centers
summary_df <- pip_df %>%
  group_by(bin_indices) %>%
  summarise(
    mean_pip = mean(pip_scores),
    eip = mean(causal_flags),
    count = n()
  ) %>%
  mutate(
    se = sqrt(eip * (1 - eip) / count),
    ci_lower = pmax(eip - 1.96 * se, 0),
    ci_upper = pmin(eip + 1.96 * se, 1),
    # Map the bin index to the corresponding bin center
    bin_center = bin_centers[as.numeric(bin_indices)]
  )


summary_df= summary_df[-which(summary_df$count<60),]
# Plot PIP Calibration with Confidence Intervals
P23 =ggplot(summary_df, aes(x =bin_center, y = eip)) +
  geom_point(color = "#1A85FF", size = 3) +  # Scatter plot of PIP vs. EIP
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.02, color = "#1A85FF") +  # CI bars
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Ideal calibration line
  labs( 
    x = "Mean PIP in Bin", y = "Empirical Inclusion Probability  ") +
  # ggtitle("fsusie  mixture per scale")+
  theme_minimal()






pip_scores <-score_susie
# 
causal_flags <-true_lab

# Define bins

bin_edges <- seq(0, 1, length.out = num_bins + 1)
bin_indices <- cut(pip_scores, breaks = bin_edges, include.lowest = TRUE, labels = FALSE)

# Calculate bin centers
bin_centers <- (bin_edges[-length(bin_edges)] + bin_edges[-1]) / 2

# Create a data frame including the bin index
pip_df <- data.frame(pip_scores, causal_flags, bin_indices)

# Compute summary statistics per bin and add bin centers
summary_df <- pip_df %>%
  group_by(bin_indices) %>%
  summarise(
    mean_pip = mean(pip_scores),
    eip = mean(causal_flags),
    count = n()
  ) %>%
  mutate(
    se = sqrt(eip * (1 - eip) / count),
    ci_lower = pmax(eip - 1.96 * se, 0),
    ci_upper = pmin(eip + 1.96 * se, 1),
    # Map the bin index to the corresponding bin center
    bin_center = bin_centers[as.numeric(bin_indices)]
  )

summary_df= summary_df[-which(summary_df$count<60),]
# Plot PIP Calibration with Confidence Intervals
P33 =ggplot(summary_df, aes(x = bin_center, y = eip)) +
  geom_point(color = "#40B0A6", size = 3) +  # Scatter plot of PIP vs. EIP
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.02, color = "#40B0A6") +  # CI bars
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Ideal calibration line
  labs( 
    x = "Mean PIP in Bin", y = "Empirical Inclusion Probability  ") +
  # ggtitle("fsusie normal mixture")+
  theme_minimal()



library(ggplot2)
library(gridExtra)
library(grid)
library(ggpubr)
library(cowplot)
library(dplyr)  

# Create text labels for the column titles
titles <- list(   textGrob(label = "Gaussian functional", gp = gpar(fontsize = 16, fontface = "bold")),
                  textGrob(label = "WGBS block", gp = gpar(fontsize = 16, fontface = "bold")),
                  textGrob(label = "WGBS distance decay", gp = gpar(fontsize = 16, fontface = "bold"))
)  

# Create text labels for the row annotations (h^2 values) using LaTeX-style expressions
h2_labels <- list(
  textGrob(expression("fSuSiE IS"), rot = 90, gp = gpar(fontsize = 14, fontface = "bold")),
  textGrob(expression("fSuSiE SPS"), rot = 90, gp = gpar(fontsize = 14, fontface = "bold")),
  textGrob(expression( "SuSiE"), rot = 90, gp = gpar(fontsize = 14, fontface = "bold"))
)


# Arrange the grid with correct alignment
P_out <- grid.arrange(
  arrangeGrob(
    textGrob(""), titles[[1]], titles[[2]], titles[[3]],  # Empty space for alignment
    ncol = 4, widths = c(0.08, 1, 1,1)
  ),
  arrangeGrob(
    h2_labels[[1]], 
    P11 + theme(legend.position = "none"), 
    P12 + theme(legend.position = "none"),
    P13 + theme(legend.position = "none"),
    
    ncol =4, widths = c(0.08, 1, 1,1)
  ),
  arrangeGrob(
    h2_labels[[2]], 
    P21 + theme(legend.position = "none"), 
    P22 + theme(legend.position = "none"),
    P23 + theme(legend.position = "none"),
    
    ncol =4, widths = c(0.08, 1, 1,1)
  ),
  arrangeGrob(
    h2_labels[[3]], 
    P31 + theme(legend.position = "none"), 
    P32 + theme(legend.position = "none"),
    P33 + theme(legend.position = "none"),
    
    ncol =4, widths = c(0.08, 1, 1,1)
  ), 
  heights = c(0.08, 1, 1, 1, 0.1)
)

# Display the final plot
P_out

ggsave(P_out , file=paste0(save_path,"pip_calibration.pdf"),
       width = 29.7,
       height = 21,
       units = "cm"
)

