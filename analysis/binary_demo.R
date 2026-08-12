library(ggplot2)
library(cowplot)
library(susieR)
library(fsusieR)
set.seed(2)

# This is a function we will use in the code chunks below to visualize
# the results of the fSuSiE fine-mapping analysis.
pip_plot <- function (fit, causal_snps = NULL, xlab = "") {
  cs_colors <- c("dodgerblue","limegreen","darkorange")
  pip <- fit$pip
  n   <- length(pip)
  dat <- data.frame(pos = 1:n,pip = pip,CS = as.character(NA))
  num_cs <- length(fit$cs)
  for (i in num_cs:1) {
    j <- fit$cs[[i]]
	dat[j,"CS"] <- paste0("CS",i)
  }
  dat_cs <- subset(dat,!is.na(CS))
  out <- ggplot(dat,aes(x = pos,y = pip)) +
    geom_point() +
	geom_point(data = dat_cs,shape = 1,size = 3,
	           mapping = aes(color = CS)) +
	scale_x_continuous(breaks = NULL) +
	scale_color_manual(values = cs_colors) + 
	labs(x = xlab,y = "PIP") +
    theme_cowplot(font_size = 12)
  if (!is.null(causal_snps)) {
	dat_causal <- dat[causal_snps,]
    out <- out + geom_point(data = dat_causal,color = "red")
  }
  return(out)
}

# Simulate a data set with 2 causal SNPs.
sigmoid <- function(x)
  1/(1 + exp(-x))
m <- 512
baseline <- rep(-2,m)
effect1 <- rep(0,m)
effect2 <- rep(0,m)
effect1[100:200] <- 0.5
effect2[400:500] <- 1
X <- N3finemapping$X
s <- apply(X,2,sd)
j <- which(s > 0.001)
X <- X[,j]
n <- nrow(X)
p <- ncol(X)
for (j in 1:p)
  X[,j] <- X[,j] - min(X[,j])
true_pos <- c(124,174)
Y <- matrix(0,n,m)
for ( i in 1:nrow(Y)){
  prob <- sigmoid(baseline +
                  X[i,true_pos[1]] * effect1 +
                  X[i,true_pos[2]] * effect2)
  Y[i,] <- rbinom(m,1,prob)
} 

# Fit an fSuSiE model
fit <- susiF(Y,X,L = 2,quantile_trans = FALSE,filter_cs = FALSE,
             maxit = 100)

# Plot the PIPs and CSs.
p1 <- pip_plot(fit,true_pos) + xlim(0,200)
print(fit$cs)
print(true_pos)

# Plot the estimated effects of the two causal SNPs.
par(mfrow = c(2,1),mar = c(4,4,2,1))
res <- get_fitted_effect(fit,l = 2,cred_band = TRUE,alpha = 0.05)
i <- which(res$cred_band["low",] > 0)
j <- which(res$cred_band["low",] > 0 & effect1 != 0)
plot(1:m,res$effect,col = "gray",pch = 20,xlab = "")
points(i,res$effect[i],col = "dodgerblue",pch = 20)
points(j,res$effect[j],col = "red",pch = 1)
res <- get_fitted_effect(fit,l = 1,cred_band = TRUE,alpha = 0.05)
i <- which(res$cred_band["low",] > 0)
j <- which(res$cred_band["low",] > 0 & effect2 != 0)
plot(1:m,res$effect,col = "gray",pch = 20,xlab = "")
points(i,res$effect[i],col = "dodgerblue",pch = 20)
points(j,res$effect[j],col = "red",pch = 1)
