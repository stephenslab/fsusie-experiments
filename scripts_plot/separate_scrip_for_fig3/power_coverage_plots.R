rm(list=ls()) 
library(dplyr) 
colors <- c("#D41159","#1A85FF","#40B0A6" )

### For PVE=10% ----
#### Calibration and power plot ----

#### Gaussian functional ---- 
path <- getwd()
load(paste( path,"/simulation/Simulation_results/check_L_accuracy_128_sd1.RData", sep=""))

df_simu <- do.call(rbind, res)

colnames(df_simu) <- c("n_cs_nps",
                       "n_effect_nps",
                       "n_false_effect_nps",
                       "mean_purity_nps",
                       "cs_size_nps",
                       "n_cs_ps",
                       "n_effect_ps",
                       "n_false_effect_ps",
                       "mean_purity_ps",
                       "cs_size_ps",
                       "Number_effect",
                       "reg_sim")
df_simu <- as.data.frame(df_simu)
df_simu <- df_simu[-which(df_simu$cs_size_ps>10),]
library(dplyr)
df_simu$power_nps <- df_simu$n_effect_nps / df_simu$Number_effect
df_simu$power_ps  <- df_simu$n_effect_ps  / df_simu$Number_effect
df_simu$t1_nps    <- 1 - df_simu$n_false_effect_nps/(df_simu$n_effect_nps + df_simu$n_false_effect_nps)
df_simu$t1_ps     <- 1 - df_simu$n_false_effect_ps/(df_simu$n_effect_ps  + df_simu$n_false_effect_ps)

mean_power_nps <- rep(NA, max(df_simu$Number_effect))
mean_power_ps  <- rep(NA, max(df_simu$Number_effect))
mean_T1_nps    <- rep(NA, max(df_simu$Number_effect))
mean_T1_ps     <- rep(NA, max(df_simu$Number_effect))

mean_cs_size_nps <- rep(NA, max(df_simu$Number_effect))
mean_cs_size_ps  <- rep(NA, max(df_simu$Number_effect))

mean_purity_nps <- rep(NA, max(df_simu$Number_effect))
mean_purity_ps  <- rep(NA, max(df_simu$Number_effect))
which_L         <- rep(NA, max(df_simu$Number_effect))

h <- 1
for (i in unique(df_simu$Number_effect)) {
  idx <- which(df_simu$Number_effect == i)
  
  mean_power_nps[h] <- mean(df_simu$power_nps[idx])
  mean_power_ps[h]  <- mean(df_simu$power_ps[idx])
  
  mean_T1_nps[h] <- mean(df_simu$t1_nps[idx])
  mean_T1_ps[h]  <- mean(df_simu$t1_ps[idx])
  
  mean_cs_size_nps[h] <- median(df_simu$cs_size_nps[idx])
  mean_cs_size_ps[h]  <- median(df_simu$cs_size_ps[idx])
  
  mean_purity_nps[h] <- mean(df_simu$mean_purity_nps[idx])
  mean_purity_ps[h]  <- mean(df_simu$mean_purity_ps[idx])
  
  which_L[h] <- i
  
  h <- h + 1
}

final_df1 <- rbind(mean_power_nps, mean_power_ps,
                   mean_T1_nps,   mean_T1_ps,
                   mean_cs_size_nps,
                   mean_cs_size_ps,
                   mean_purity_nps,
                   mean_purity_ps,
                   which_L)

t_names <- rownames(final_df1)
final_df1 <- data.frame(t(final_df1))
colnames(final_df1) <- t_names

df_plot <- data.frame(power      = c(final_df1$mean_power_ps, final_df1$mean_power_nps),
                      T1_error   = c(final_df1$mean_T1_ps,   final_df1$mean_T1_nps),
                      cs_size    = c(final_df1$mean_cs_size_ps, final_df1$mean_cs_size_nps),
                      mean_purity= c(final_df1$mean_purity_ps, final_df1$mean_purity_nps),
                      prior      = as.factor(c(rep("SPSP", nrow(final_df1)),
                                               rep("ISP",  nrow(final_df1)))),
                      L          = factor(rep(final_df1$which_L, 2)))
df_plot <- df_plot[complete.cases(df_plot),]

## --- Replace normal approximation CI with exact CI from binom.test ---
# Define functions that take a (mean) proportion and an effective replicate number
# and return the lower or upper bound of the exact binom.test confidence interval.
exact_ci_low <- function(p, n) {
  if(is.na(p)){
    return(NA)
  }
  successes <- round(p * n)
  ci <- binom.test(successes, n)$conf.int
  return(ci[1])
}
exact_ci_up <- function(p, n) {
  if(is.na(p)){
    return(NA)
  }
  successes <- round(p * n)
  ci <- binom.test(successes, n)$conf.int
  return(ci[2])
}

# Here, we use n_rep as the effective number of simulation replicates for each group.
# (In your original code, n_rep was defined via rep(c(table(df_simu$Number_effect)), each=2).)
n_rep <- rep(as.numeric(table(df_simu$Number_effect)), each=2)

# Now compute the error bar endpoints for both power and T1_error using binom.test.
df_plot$pw_er_low <- mapply(exact_ci_low, df_plot$power, n_rep)
df_plot$pw_er_up  <- mapply(exact_ci_up,  df_plot$power, n_rep)
df_plot$t1_er_low <- mapply(exact_ci_low, df_plot$T1_error, n_rep)
df_plot$t1_er_up  <- mapply(exact_ci_up,  df_plot$T1_error, n_rep)

### (The rest of your code for plotting remains the same.)
path <- getwd()
load(paste0(path, "/simulation/Simulation_script/script/additional_simualtion_for_fig3_panel_D_susie_pc_calibration_power/check_L_susie_128_sd1.RData"))

df_susie <- data.frame(do.call(rbind, res))
colnames(df_susie) <- c("Number_effect", "n_cs", "n_effect", "n_false_effect",
                        "purity",
                        "cs_size")
df_simu <- df_susie
df_simu$power <- df_simu$n_effect / df_simu$Number_effect 
df_simu$t1    <- 1 - df_simu$n_false_effect / (df_simu$n_effect + df_simu$n_false_effect)

mean_power_susie <- rep(NA, max(df_simu$Number_effect))
mean_T1_susie    <- rep(NA, max(df_simu$Number_effect))

mean_purity    <- rep(NA, max(df_simu$Number_effect))
mean_cs    <- rep(NA, max(df_simu$Number_effect))
which_L         <- rep(NA, max(df_simu$Number_effect))

h <- 1
for (i in unique(df_simu$Number_effect)) {
  idx <- which(df_simu$Number_effect == i)
  
  mean_power_susie[h] <- mean(df_simu$power[idx])
  mean_T1_susie[h]    <- mean(df_simu$t1[idx])
  mean_purity [h]     <- mean(df_simu$purity[idx])
  mean_cs [h]         <- mean(df_simu$cs_size[idx])
  which_L[h] <- i
  h <- h + 1
}

df <- data.frame(power    = mean_power_susie,
                 T1_error = mean_T1_susie,
                 L        = unique(df_simu$Number_effect),
                 mean_purity   =  mean_purity,
                 cs_size  = mean_cs)

plot(df$L, df$power)
plot(df$L, df$T1_error)

# For the SuSiE data, we again replace the normal approximation error bars.
n_rep <- as.numeric(table(df_simu$Number_effect))  # here one value per group
df$pw_er_low <- mapply(exact_ci_low, df$power, n_rep)
df$pw_er_up  <- mapply(exact_ci_up,  df$power, n_rep)
df$t1_er_low <- mapply(exact_ci_low, df$T1_error, n_rep)
df$t1_er_up  <- mapply(exact_ci_up,  df$T1_error, n_rep)
df$L <- as.factor(df$L)
df$prior <- rep("SuSiE", nrow(df))



dfp <- rbind(df_plot[, c("power", "T1_error", "prior", "L",
                         "pw_er_up", "pw_er_low", "t1_er_up", "t1_er_low",
                         "cs_size","mean_purity")],
             df[, c("power", "T1_error", "prior", "L",
                    "pw_er_up", "pw_er_low", "t1_er_up", "t1_er_low",
                    
                    "cs_size","mean_purity")])
library(ggplot2)

P11 <- ggplot(dfp, aes(x = L, y = power, col = prior)) +
  geom_point(position = position_dodge(.9), size = 2) +
  geom_errorbar(aes(ymin = pw_er_low, ymax = pw_er_up), width = .2,
                position = position_dodge(.9)) +
  ylim(c(0, 1)) +
  ylab("Power") +
  theme_linedraw() +
  scale_color_manual(values = colors) +
  scale_x_discrete(labels = c(1, "", 3, "", 5, "", 7, "", 9, "", 11, "", 13, "", 15, "", 17, "", 19, ""))

P11_t1 <- ggplot(dfp, aes(x = L, y = T1_error, col = prior)) +
  geom_hline(yintercept = 0.95, color = "black", linewidth = 2) +
  geom_point(position = position_dodge(.9), size = 2) +
  geom_errorbar(aes(ymin = t1_er_low, ymax = t1_er_up), width = .2,
                position = position_dodge(.9)) +
  ylim(c(0.7, 1.04)) +
  ylab("Coverage") +
  theme_linedraw() +
  scale_color_manual(values = colors) +
  scale_x_discrete(labels = c(1, "", 3, "", 5, "", 7, "", 9, "", 11, "", 13, "", 15, "", 17, "", 19, ""))

P11_t1
#################################
## (2) Block Section          ##
#################################
path <- getwd()
load(paste(path, "/simulation/Simulation_results/block_L_accuracy_128_sd1.RData", sep=""))
df_simu <- do.call(rbind, res)
colnames(df_simu) <- c("n_cs_nps",
                       "n_effect_nps",
                       "n_false_effect_nps",
                       "mean_purity_nps",
                       "cs_size_nps",
                       "n_cs_ps",
                       "n_effect_ps",
                       "n_false_effect_ps",
                       "mean_purity_ps",
                       "cs_size_ps",
                       "Number_effect",
                       "reg_sim")
df_simu <- as.data.frame(df_simu)
df_simu <- df_simu[-which(df_simu$cs_size_ps > 10), ]
df_simu$power_nps <- df_simu$n_effect_nps / df_simu$Number_effect
df_simu$power_ps  <- df_simu$n_effect_ps  / df_simu$Number_effect
df_simu$t1_nps    <- 1 - df_simu$n_false_effect_nps / (df_simu$n_effect_nps + df_simu$n_false_effect_nps)
df_simu$t1_ps     <- 1 - df_simu$n_false_effect_ps  / (df_simu$n_effect_ps  + df_simu$n_false_effect_ps)

# Aggregate block simulation results by Number_effect
mean_power_nps <- rep(NA, max(df_simu$Number_effect))
mean_power_ps  <- rep(NA, max(df_simu$Number_effect))
mean_T1_nps    <- rep(NA, max(df_simu$Number_effect))
mean_T1_ps     <- rep(NA, max(df_simu$Number_effect))
mean_cs_size_nps <- rep(NA, max(df_simu$Number_effect))
mean_cs_size_ps  <- rep(NA, max(df_simu$Number_effect))
mean_purity_nps <- rep(NA, max(df_simu$Number_effect))
mean_purity_ps  <- rep(NA, max(df_simu$Number_effect))
which_L         <- rep(NA, max(df_simu$Number_effect))
h <- 1
for (i in unique(df_simu$Number_effect)) {
  idx <- which(df_simu$Number_effect == i)
  mean_power_nps[h] <- mean(df_simu$power_nps[idx])
  mean_power_ps[h]  <- mean(df_simu$power_ps[idx])
  mean_T1_nps[h]    <- mean(df_simu$t1_nps[idx])
  mean_T1_ps[h]     <- mean(df_simu$t1_ps[idx])
  mean_cs_size_nps[h] <- median(df_simu$cs_size_nps[idx])
  mean_cs_size_ps[h]  <- median(df_simu$cs_size_ps[idx])
  mean_purity_nps[h]  <- mean(df_simu$mean_purity_nps[idx])
  mean_purity_ps[h]   <- mean(df_simu$mean_purity_ps[idx])
  which_L[h] <- i
  h <- h + 1
}
final_df1 <- rbind(mean_power_nps, mean_power_ps,
                   mean_T1_nps,   mean_T1_ps,
                   mean_cs_size_nps,
                   mean_cs_size_ps,
                   mean_purity_nps,
                   mean_purity_ps,
                   which_L)
t_names <- rownames(final_df1)
final_df1 <- data.frame(t(final_df1))
colnames(final_df1) <- t_names

df_plot <- data.frame(power       = c(final_df1$mean_power_ps, final_df1$mean_power_nps),
                      T1_error    = c(final_df1$mean_T1_ps,   final_df1$mean_T1_nps),
                      cs_size     = c(final_df1$mean_cs_size_ps, final_df1$mean_cs_size_nps),
                      mean_purity = c(final_df1$mean_purity_ps,  final_df1$mean_purity_nps),
                      prior       = as.factor(c(rep("SPSP", nrow(final_df1)),
                                                rep("ISP",  nrow(final_df1)))),
                      L           = factor(rep(final_df1$which_L, 2)))
df_plot <- df_plot[complete.cases(df_plot), ]
# Use exact binom.test CI for block error bars
n_rep <- rep(as.numeric(table(df_simu$Number_effect)), each = 2)
df_plot$pw_er_low <- mapply(exact_ci_low, df_plot$power, n_rep)
df_plot$pw_er_up  <- mapply(exact_ci_up,  df_plot$power, n_rep)
df_plot$t1_er_low <- mapply(exact_ci_low, df_plot$T1_error, n_rep)
df_plot$t1_er_up  <- mapply(exact_ci_up,  df_plot$T1_error, n_rep)

# Process SuSiE block results
path <- getwd()
load(paste0(path, "/simulation/Simulation_script/script/additional_simualtion_for_fig3_panel_D_susie_pc_calibration_power/check_L_susie_128_block_sd1.RData"))

df_susie <- data.frame(do.call(rbind, res))
colnames(df_susie) <- c("Number_effect", "n_cs", "n_effect", "n_false_effect",
                        "purity",
                        "cs_size")
df_simu <- df_susie
df_simu$power <- df_simu$n_effect / df_simu$Number_effect 
df_simu$t1    <- 1 - df_simu$n_false_effect / (df_simu$n_effect + df_simu$n_false_effect)

mean_power_susie <- rep(NA, max(df_simu$Number_effect))
mean_T1_susie    <- rep(NA, max(df_simu$Number_effect))

mean_purity    <- rep(NA, max(df_simu$Number_effect))
mean_cs    <- rep(NA, max(df_simu$Number_effect))
which_L         <- rep(NA, max(df_simu$Number_effect))

h <- 1
for (i in unique(df_simu$Number_effect)) {
  idx <- which(df_simu$Number_effect == i)
  
  mean_power_susie[h] <- mean(df_simu$power[idx])
  mean_T1_susie[h]    <- mean(df_simu$t1[idx])
  mean_purity [h]     <- mean(df_simu$purity[idx])
  mean_cs [h]         <- mean(df_simu$cs_size[idx])
  which_L[h] <- i
  h <- h + 1
}

df <- data.frame(power    = mean_power_susie,
                 T1_error = mean_T1_susie,
                 L        = unique(df_simu$Number_effect),
                 mean_purity  =  mean_purity,
                 cs_size  = mean_cs)
# For SuSiE block, use exact CI functions
n_rep <- as.numeric(table(df_simu$Number_effect))
df$pw_er_low <- mapply(exact_ci_low, df$power, n_rep)
df$pw_er_up  <- mapply(exact_ci_up,  df$power, n_rep)
df$t1_er_low <- mapply(exact_ci_low, df$T1_error, n_rep)
df$t1_er_up  <- mapply(exact_ci_up,  df$T1_error, n_rep)
df$L <- as.factor(df$L)
df$prior <- rep("SuSiE", nrow(df))
dfp_block <-  rbind(df_plot[, c("power", "T1_error", "prior", "L",
                                "pw_er_up", "pw_er_low", "t1_er_up", "t1_er_low",
                                "cs_size","mean_purity")],
                    df[, c("power", "T1_error", "prior", "L",
                           "pw_er_up", "pw_er_low", "t1_er_up", "t1_er_low",
                           
                           "cs_size","mean_purity")])

# Plot block section results
library(ggplot2)
P21 <- ggplot(dfp_block, aes(x = L, y = power, col = prior)) +
  geom_point(position = position_dodge(.9), size = 2) +
  geom_errorbar(aes(ymin = pw_er_low, ymax = pw_er_up), width = .2,
                position = position_dodge(.9)) +
  ylim(c(0, 1)) +
  ylab("Power") +
  theme_linedraw() +
  scale_color_manual(values = colors) +
  scale_x_discrete(labels = c(1, "", 3, "", 5, "", 7, "", 9, "", 11, "", 13, "", 15, "", 17, "", 19, ""))
P21

P21_t1 <- ggplot(dfp_block, aes(x = L, y = T1_error, col = prior)) +
  geom_hline(yintercept = 0.95, color = "black", linewidth = 2) +
  geom_point(position = position_dodge(.9), size = 2) +
  geom_errorbar(aes(ymin = t1_er_low, ymax = t1_er_up), width = .2,
                position = position_dodge(.9)) +
  ylim(c(0.7, 1.04)) +
  ylab("Coverage") +
  theme_linedraw() +
  scale_color_manual(values = colors) +
  scale_x_discrete(labels = c(1, "", 3, "", 5, "", 7, "", 9, "", 11, "", 13, "", 15, "", 17, "", 19, ""))
P21_t1

#################################
## (3) Distance Decay Section ##
#################################
path <- getwd()
load(paste(path, "/simulation/Simulation_results/dist_decay_L_accuracy_128_sd1.RData", sep=""))
df_simu <- do.call(rbind, res)
colnames(df_simu) <- c("n_cs_nps",
                       "n_effect_nps",
                       "n_false_effect_nps",
                       "mean_purity_nps",
                       "cs_size_nps",
                       "n_cs_ps",
                       "n_effect_ps",
                       "n_false_effect_ps",
                       "mean_purity_ps",
                       "cs_size_ps",
                       "Number_effect",
                       "reg_sim")
df_simu <- as.data.frame(df_simu)
df_simu <- df_simu[-which(df_simu$cs_size_ps > 10), ]
df_simu$power_nps <- df_simu$n_effect_nps / df_simu$Number_effect
df_simu$power_ps  <- df_simu$n_effect_ps  / df_simu$Number_effect
df_simu$t1_nps    <- 1 - df_simu$n_false_effect_nps / (df_simu$n_effect_nps + df_simu$n_false_effect_nps)
df_simu$t1_ps     <- 1 - df_simu$n_false_effect_ps  / (df_simu$n_effect_ps  + df_simu$n_false_effect_ps)

# Aggregate distance decay simulation results by Number_effect
mean_power_nps <- rep(NA, max(df_simu$Number_effect))
mean_power_ps  <- rep(NA, max(df_simu$Number_effect))
mean_T1_nps    <- rep(NA, max(df_simu$Number_effect))
mean_T1_ps     <- rep(NA, max(df_simu$Number_effect))
mean_cs_size_nps <- rep(NA, max(df_simu$Number_effect))
mean_cs_size_ps  <- rep(NA, max(df_simu$Number_effect))
mean_purity_nps <- rep(NA, max(df_simu$Number_effect))
mean_purity_ps  <- rep(NA, max(df_simu$Number_effect))
which_L         <- rep(NA, max(df_simu$Number_effect))
h <- 1
for (i in unique(df_simu$Number_effect)) {
  idx <- which(df_simu$Number_effect == i)
  mean_power_nps[h] <- mean(df_simu$power_nps[idx])
  mean_power_ps[h]  <- mean(df_simu$power_ps[idx])
  mean_T1_nps[h]    <- mean(df_simu$t1_nps[idx])
  mean_T1_ps[h]     <- mean(df_simu$t1_ps[idx])
  mean_cs_size_nps[h] <- median(df_simu$cs_size_nps[idx])
  mean_cs_size_ps[h]  <- median(df_simu$cs_size_ps[idx])
  mean_purity_nps[h]  <- mean(df_simu$mean_purity_nps[idx])
  mean_purity_ps[h]   <- mean(df_simu$mean_purity_ps[idx])
  which_L[h] <- i
  h <- h + 1
}
final_df1 <- rbind(mean_power_nps, mean_power_ps,
                   mean_T1_nps,   mean_T1_ps,
                   mean_cs_size_nps,
                   mean_cs_size_ps,
                   mean_purity_nps,
                   mean_purity_ps,
                   which_L)
t_names <- rownames(final_df1)
final_df1 <- data.frame(t(final_df1))
colnames(final_df1) <- t_names

df_plot <- data.frame(power       = c(final_df1$mean_power_ps, final_df1$mean_power_nps),
                      T1_error    = c(final_df1$mean_T1_ps,   final_df1$mean_T1_nps),
                      cs_size     = c(final_df1$mean_cs_size_ps, final_df1$mean_cs_size_nps),
                      mean_purity = c(final_df1$mean_purity_ps,  final_df1$mean_purity_nps),
                      prior       = as.factor(c(rep("SPSP", nrow(final_df1)),
                                                rep("ISP",  nrow(final_df1)))),
                      L           = factor(rep(final_df1$which_L, 2)))
df_plot <- df_plot[complete.cases(df_plot), ]
# Use exact CI for distance decay error bars
n_rep <- rep(as.numeric(table(df_simu$Number_effect)), each = 2)
df_plot$pw_er_low <- mapply(exact_ci_low, df_plot$power, n_rep)
df_plot$pw_er_up  <- mapply(exact_ci_up,  df_plot$power, n_rep)
df_plot$t1_er_low <- mapply(exact_ci_low, df_plot$T1_error, n_rep)
df_plot$t1_er_up  <- mapply(exact_ci_up,  df_plot$T1_error, n_rep)

# Process SuSiE distance decay results
path <- getwd()
load(paste0(path, "/simulation/Simulation_script/script/additional_simualtion_for_fig3_panel_D_susie_pc_calibration_power/check_L_susie_128_distdecay_sd1.RData"))

df_susie <- data.frame(do.call(rbind, res))
colnames(df_susie) <- c("Number_effect", "n_cs", "n_effect", "n_false_effect",
                        "purity",
                        "cs_size")
df_simu <- df_susie
df_simu$power <- df_simu$n_effect / df_simu$Number_effect 
df_simu$t1    <- 1 - df_simu$n_false_effect / (df_simu$n_effect + df_simu$n_false_effect)

mean_power_susie <- rep(NA, max(df_simu$Number_effect))
mean_T1_susie    <- rep(NA, max(df_simu$Number_effect))

mean_purity    <- rep(NA, max(df_simu$Number_effect))
mean_cs    <- rep(NA, max(df_simu$Number_effect))
which_L         <- rep(NA, max(df_simu$Number_effect))

h <- 1
for (i in unique(df_simu$Number_effect)) {
  idx <- which(df_simu$Number_effect == i)
  
  mean_power_susie[h] <- mean(df_simu$power[idx])
  mean_T1_susie[h]    <- mean(df_simu$t1[idx])
  mean_purity [h]     <- mean(df_simu$purity[idx])
  mean_cs [h]         <- mean(df_simu$cs_size[idx])
  which_L[h] <- i
  h <- h + 1
}

df <- data.frame(power    = mean_power_susie,
                 T1_error = mean_T1_susie,
                 L        = unique(df_simu$Number_effect),
                 mean_purity   =  mean_purity,
                 cs_size  = mean_cs)
# For SuSiE distance decay, use exact CI functions
n_rep <- as.numeric(table(df_simu$Number_effect))
df$pw_er_low <- mapply(exact_ci_low, df$power, n_rep)
df$pw_er_up  <- mapply(exact_ci_up,  df$power, n_rep)
df$t1_er_low <- mapply(exact_ci_low, df$T1_error, n_rep)
df$t1_er_up  <- mapply(exact_ci_up,  df$T1_error, n_rep)
df$L <- as.factor(df$L)
df$prior <- rep("SuSiE", nrow(df))
dfp_decay <- rbind(df_plot[, c("power", "T1_error", "prior", "L",
                               "pw_er_up", "pw_er_low", "t1_er_up", "t1_er_low",
                               "cs_size","mean_purity")],
                   df[, c("power", "T1_error", "prior", "L",
                          "pw_er_up", "pw_er_low", "t1_er_up", "t1_er_low",
                          
                          "cs_size","mean_purity")])
# Plot distance decay section results
P31 <- ggplot(dfp_decay, aes(x = L, y = power, col = prior)) +
  geom_point(position = position_dodge(.9), size = 2) +
  geom_errorbar(aes(ymin = pw_er_low, ymax = pw_er_up), width = .2,
                position = position_dodge(.9)) +
  ylim(c(0, 1)) +
  ylab("Power") +
  theme_linedraw() +
  scale_color_manual(values = colors) +
  scale_x_discrete(labels = c(1, "", 3, "", 5, "", 7, "", 9, "", 11, "", 13, "", 15, "", 17, "", 19, ""))
P31

P31_t1 <- ggplot(dfp_decay, aes(x = L, y = T1_error, col = prior)) +
  geom_hline(yintercept = 0.95, color = "black", linewidth = 2) +
  geom_point(position = position_dodge(.9), size = 2) +
  geom_errorbar(aes(ymin = t1_er_low, ymax = t1_er_up), width = .2,
                position = position_dodge(.9)) +
  ylim(c(0.7, 1.04)) +
  ylab("Coverage") +
  theme_linedraw() +
  scale_color_manual(values = colors) +
  scale_x_discrete(labels = c(1, "", 3, "", 5, "", 7, "", 9, "", 11, "", 13, "", 15, "", 17, "", 19, ""))
P31_t1

#################################
## Final Plot Arrangement     ##
#################################
library(gridExtra)
library(cowplot)
library(cowplot)
library(dplyr)
library(susieR)
library(ggpubr)
library(gridExtra)
library(grid)
titles <- list(textGrob(label = "Gaussian functional", gp = gpar(fontsize = 16, fontface = "bold")),
               textGrob(label = "WGBS block",       gp = gpar(fontsize = 16, fontface = "bold")),
               textGrob(label = "WGBS distance decay", gp = gpar(fontsize = 16, fontface = "bold")))
legend_plot <- get_legend(P31) + 
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.title = element_blank())

# Example: arranging Gaussian functional plots (P11 & P11_t1)
P_out <- grid.arrange(
  arrangeGrob(
    titles[[1]], titles[[2]], titles[[3]],
    ncol = 3, widths = c(1, 1, 1)
  ),
  arrangeGrob(
    P11 + theme(legend.position = "none"),
    P21 + theme(legend.position = "none"),
    P31 + theme(legend.position = "none"),
    ncol = 3, widths = c(1, 1, 1)
  ),
  heights = c(0.08, 1)
)
P_out

# Save combined power plot
save_path <- paste0(getwd(), "/plot/fig3_separate_panel/")
ggsave(P_out, file = paste0(save_path, "power.pdf"),
       width = 29.7, height = 10.5, units = "cm")

# Arrange and save coverage plot similarly
P_out_cov <- grid.arrange(
  arrangeGrob(
    titles[[1]], titles[[2]], titles[[3]],
    ncol = 3, widths = c(1, 1, 1)
  ),
  arrangeGrob(
    P11_t1 + theme(legend.position = "none"),
    P21_t1 + theme(legend.position = "none"),
    P31_t1 + theme(legend.position = "none"),
    ncol = 3, widths = c(1, 1, 1)
  ),
  heights = c(0.08, 1)
)
ggsave(P_out_cov, file = paste0(save_path, "coverage.pdf"),
       width = 29.7, height = 10.5, units = "cm")









library(dplyr)
decay = dfp_decay[which(as.numeric(dfp_decay$L)<15 ),] %>%
  group_by_("prior")%>% 
  summarise(power       = mean(power,na.rm=TRUE),
            coverage    = mean(T1_error, na.rm=TRUE),
            purity      = mean(mean_purity, na.rm=TRUE),
            purity_up   = mean(mean_purity, na.rm=TRUE)+ mean( 1.96*sqrt( var(mean_purity,na.rm=TRUE))),# /300)),
            purity_low  = mean(mean_purity, na.rm=TRUE)- mean( 1.96*sqrt( var(mean_purity,na.rm=TRUE))),# /300)),
            cs_size_est     = mean(cs_size,na.rm=TRUE),
            cs_size_up  =cs_size_est +1.96*   mean( 1.96*sqrt( var(cs_size,na.rm=TRUE))),# /300)),
            cs_size_low =cs_size_est -1.96*    mean( 1.96*sqrt( var(cs_size,na.rm=TRUE)))# /300)),
            )
decay 
block= dfp_block[which(as.numeric(dfp_block$L)<15),] %>%
  group_by_("prior")%>% 
  summarise(power       = mean(power,na.rm=TRUE),
            coverage    = mean(T1_error, na.rm=TRUE),
            purity      = mean(mean_purity, na.rm=TRUE),
            purity_up   = mean(mean_purity, na.rm=TRUE)+ mean( 1.96*sqrt( var(mean_purity,na.rm=TRUE))),# /300)),
            purity_low  = mean(mean_purity, na.rm=TRUE)- mean( 1.96*sqrt( var(mean_purity,na.rm=TRUE))),# /300)),
            cs_size_est     = mean(cs_size,na.rm=TRUE),
            cs_size_up  =cs_size_est +1.96*   mean( 1.96*sqrt( var(cs_size,na.rm=TRUE))),# /300)),
            cs_size_low =cs_size_est -1.96*    mean( 1.96*sqrt( var(cs_size,na.rm=TRUE)))# /300)),
  )

gaussian= dfp[which(as.numeric(dfp$L)<15 ),] %>%
  group_by_("prior")%>% 
  summarise(power       = mean(power,na.rm=TRUE),
            coverage    = mean(T1_error, na.rm=TRUE),
            purity      = mean(mean_purity, na.rm=TRUE),
            purity_up   = mean(mean_purity, na.rm=TRUE)+ mean( 1.96*sqrt( var(mean_purity,na.rm=TRUE))),# /300)),
            purity_low  = mean(mean_purity, na.rm=TRUE)- mean( 1.96*sqrt( var(mean_purity,na.rm=TRUE))),# /300)),
            cs_size_est     = mean(cs_size,na.rm=TRUE),
            cs_size_up  =cs_size_est +1.96*   mean( 1.96*sqrt( var(cs_size,na.rm=TRUE))),# /300)),
            cs_size_low =cs_size_est -1.96*    mean( 1.96*sqrt( var(cs_size,na.rm=TRUE)))# /300)),
  )


gaussian$sim=  "Gaussian"

block$sim = " WGBS block"
decay $sim = " WGBS decay"

df_panel_C= rbind(gaussian,
                  block,
                  decay)

#### Neww  ----
df_panel_C$sim <- factor(df_panel_C$sim, levels = c("Gaussian", " WGBS block", " WGBS decay"))

P_power = ggplot( df_panel_C,
        aes(y=power, x=sim, col=prior))+
  geom_point(position = position_dodge(.5), size = 2) +
  theme_linedraw() +
   theme(legend.position = "none")+
  ylab("Power")+
  xlab("")+
  
  
  scale_color_manual(values=c(  "gold","magenta","dodgerblue" ))

P_power
ggsave(P_power, file = paste0(save_path, "power_1_15.pdf"),
       width = 10.5, height = 10.5, units = "cm")
P_coverage =  ggplot( df_panel_C,
        aes(y=coverage, x=sim, col=prior))+
  geom_point(position = position_dodge(.5), size = 2) +
  theme_linedraw() +
  theme(legend.position = "none")+
  ylab("Coverage")+
  xlab("")+
  
  scale_color_manual(values=c(  "gold","magenta","dodgerblue" ))
P_coverage
ggsave(P_coverage, file = paste0(save_path, "coverage_1_15.pdf"),
       width = 10.5, height = 10.5, units = "cm")
  
P_cs_size =  ggplot( df_panel_C,
                      aes(y=cs_size_est, x=sim, col=prior))+
  geom_point(position = position_dodge(.5), size = 2) +
  
  ylab("CS size")+
  xlab("")+
  theme_linedraw() +
  theme(legend.position = "none")+
  geom_errorbar(aes(ymin = cs_size_low, ymax =cs_size_up), width = .2,
                position = position_dodge(.5)) +
  scale_color_manual(values=c(  "gold","magenta","dodgerblue" ))
P_cs_size
ggsave(P_cs_size, file = paste0(save_path, "cs_size_1_15.pdf"),
       width = 10.5, height = 10.5, units = "cm")
P_purity =  ggplot( df_panel_C,
                     aes(y=purity, x=sim, col=prior))+
  geom_point(position = position_dodge(.5), size = 2) +
  theme(legend.position = "none")+
  ylab("Purity")+
  xlab("")+
  theme_linedraw() +
  theme(legend.position = "none")+
  geom_errorbar(aes(ymin = purity_low, ymax =purity_up), width = .2,
                position = position_dodge(.5)) +
  scale_color_manual(values=c(  "gold","magenta","dodgerblue" ))
P_purity

ggsave(P_purity, file = paste0(save_path, "purity_1_15.pdf"),
       width = 10.5, height = 10.5, units = "cm")

