# Load necessary library
library(ggplot2)
library(dplyr)


path <- getwd()
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


summary_df= summary_df[-which(summary_df$count<30),]
# Plot PIP Calibration with Confidence Intervals
P11 =ggplot(summary_df, aes(x = mean_pip, y = eip)) +
  geom_point(color = "blue", size = 3) +  # Scatter plot of PIP vs. EIP
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.02, color = "blue") +  # CI bars
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Ideal calibration line
  labs(title = "PIP Calibration Plot with Confidence Intervals",
       x = "Mean PIP in Bin", y = "Empirical Inclusion Probability (EIP)") +
  ggtitle("fsusie normal mixture")+
  theme_minimal()




pip_scores <- score_sp_fsusie
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


summary_df= summary_df[-which(summary_df$count<30),]
# Plot PIP Calibration with Confidence Intervals
P21 =ggplot(summary_df, aes(x = mean_pip, y = eip)) +
  geom_point(color = "blue", size = 3) +  # Scatter plot of PIP vs. EIP
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.02, color = "blue") +  # CI bars
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Ideal calibration line
  labs(title = "PIP Calibration Plot with Confidence Intervals",
       x = "Mean PIP in Bin", y = "Empirical Inclusion Probability (EIP)") +
  ggtitle("fsusie  mixture per scale")+
  theme_minimal()






pip_scores <-score_susie
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


summary_df= summary_df[-which(summary_df$count<40),]
# Plot PIP Calibration with Confidence Intervals
P31 =ggplot(summary_df, aes(x = mean_pip, y = eip)) +
  geom_point(color = "blue", size = 3) +  # Scatter plot of PIP vs. EIP
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.02, color = "blue") +  # CI bars
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Ideal calibration line
  labs(title = "PIP Calibration Plot with Confidence Intervals",
       x = "Mean PIP in Bin", y = "Empirical Inclusion Probability (EIP)") +
  ggtitle("fsusie normal mixture")+
  theme_minimal()



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


summary_df= summary_df[-which(summary_df$count<30),]
# Plot PIP Calibration with Confidence Intervals
P12 =ggplot(summary_df, aes(x = mean_pip, y = eip)) +
  geom_point(color = "blue", size = 3) +  # Scatter plot of PIP vs. EIP
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.02, color = "blue") +  # CI bars
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Ideal calibration line
  labs(title = "PIP Calibration Plot with Confidence Intervals",
       x = "Mean PIP in Bin", y = "Empirical Inclusion Probability (EIP)") +
  ggtitle("fsusie normal mixture")+
  theme_minimal()

  
  

pip_scores <- score_sp_fsusie
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


summary_df= summary_df[-which(summary_df$count<30),]
# Plot PIP Calibration with Confidence Intervals
P22 =ggplot(summary_df, aes(x = mean_pip, y = eip)) +
  geom_point(color = "blue", size = 3) +  # Scatter plot of PIP vs. EIP
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.02, color = "blue") +  # CI bars
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Ideal calibration line
  labs(title = "PIP Calibration Plot with Confidence Intervals",
       x = "Mean PIP in Bin", y = "Empirical Inclusion Probability (EIP)") +
  ggtitle("fsusie  mixture per scale")+
  theme_minimal()

  
    
  
  
  
  pip_scores <-score_susie
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
  
  
  summary_df= summary_df[-which(summary_df$count<30),]
  # Plot PIP Calibration with Confidence Intervals
  P32 =ggplot(summary_df, aes(x = mean_pip, y = eip)) +
    geom_point(color = "blue", size = 3) +  # Scatter plot of PIP vs. EIP
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.02, color = "blue") +  # CI bars
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Ideal calibration line
    labs(title = "PIP Calibration Plot with Confidence Intervals",
         x = "Mean PIP in Bin", y = "Empirical Inclusion Probability (EIP)") +
    ggtitle("fsusie normal mixture")+
  theme_minimal()
  
  
  
  
  
  
  
  
  
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
  
  
  summary_df= summary_df[-which(summary_df$count<30),]
  # Plot PIP Calibration with Confidence Intervals
  P13 =ggplot(summary_df, aes(x = mean_pip, y = eip)) +
    geom_point(color = "blue", size = 3) +  # Scatter plot of PIP vs. EIP
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.02, color = "blue") +  # CI bars
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Ideal calibration line
    labs(title = "PIP Calibration Plot with Confidence Intervals",
         x = "Mean PIP in Bin", y = "Empirical Inclusion Probability (EIP)") +
    ggtitle("fsusie normal mixture")+
    theme_minimal()
  
  
  
  
  pip_scores <- score_sp_fsusie
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
  
  
  summary_df= summary_df[-which(summary_df$count<30),]
  # Plot PIP Calibration with Confidence Intervals
  P23 =ggplot(summary_df, aes(x = mean_pip, y = eip)) +
    geom_point(color = "blue", size = 3) +  # Scatter plot of PIP vs. EIP
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.02, color = "blue") +  # CI bars
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Ideal calibration line
    labs(title = "PIP Calibration Plot with Confidence Intervals",
         x = "Mean PIP in Bin", y = "Empirical Inclusion Probability (EIP)") +
    ggtitle("fsusie  mixture per scale")+
    theme_minimal()
  
  
  
  
  
  
  pip_scores <-score_susie
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
  
  
  summary_df= summary_df[-which(summary_df$count<30),]
  # Plot PIP Calibration with Confidence Intervals
  P33 =ggplot(summary_df, aes(x = mean_pip, y = eip)) +
    geom_point(color = "blue", size = 3) +  # Scatter plot of PIP vs. EIP
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.02, color = "blue") +  # CI bars
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Ideal calibration line
    labs(title = "PIP Calibration Plot with Confidence Intervals",
         x = "Mean PIP in Bin", y = "Empirical Inclusion Probability (EIP)") +
    ggtitle("fsusie normal mixture")+
    theme_minimal()

  
  
  
  
  
  
   
   