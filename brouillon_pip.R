# Load necessary library
library(ggplot2)

# Simulate data: PIP scores and ground-truth causal indicators
set.seed(123)
n <- 1000
pip_scores <- score_fsusie#score_susie #
# 
causal_flags <-true_lab

# Define bins
num_bins <- 10
bin_edges <- seq(0, 1, length.out = num_bins + 1)
bin_indices <- cut(pip_scores, breaks = bin_edges, include.lowest = TRUE, labels = FALSE)

# Compute Mean PIP and Empirical Inclusion Probability (EIP) per bin
pip_df <- data.frame(pip_scores, causal_flags, bin_indices)

# Compute summary statistics per bin
summary_df <- pip_df %>%
  group_by(bin_indices) %>%
  summarise(
    mean_pip = mean(pip_scores),
    eip = mean(causal_flags),
    count = n()
  ) %>%
  mutate(
    se = sqrt(eip * (1 - eip) / count),  # Standard error for proportion
    ci_lower = pmax(eip - 1.96 * se, 0), # 95% CI lower bound
    ci_upper = pmin(eip + 1.96 * se, 1)  # 95% CI upper bound
  )

# Plot PIP Calibration with Confidence Intervals
ggplot(summary_df, aes(x = mean_pip, y = eip)) +
  geom_point(color = "blue", size = 3) +  # Scatter plot of PIP vs. EIP
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.02, color = "blue") +  # CI bars
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Ideal calibration line
  labs(title = "PIP Calibration Plot with Confidence Intervals",
       x = "Mean PIP in Bin", y = "Empirical Inclusion Probability (EIP)") +
  theme_minimal()
