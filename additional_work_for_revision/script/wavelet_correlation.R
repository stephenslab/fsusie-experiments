## ------------------------------------------------------------
## Gaussian process -> wavelet transform -> correlation
## Uses wavethresh
## ------------------------------------------------------------
path= getwd ()
library(wavethresh)
library(MASS)

set.seed(1)

## -------------------------
## 1. Simulate GP samples
## -------------------------
n  <- 2000        # number of GP draws (rows)
S  <- 10  # grid length (must be power of 2)
ell <- 200 
sigma2 <- 1

x <- 1: 2^S
D <- as.matrix(dist(x))

# Exponential GP covariance
Sigma <- sigma2 * exp(-D / ell)

# Each row is one GP draw
Y <- mvrnorm(n = n, mu = rep(0,  2^S), Sigma = Sigma)

## -------------------------
## 2. Correlation: original domain
## -------------------------

# Correlation between locations (empirical GP correlation)
cor_time <- cor(Y)

orig_corr_mean <- mean(abs(cor_time[upper.tri(cor_time)]))

## -------------------------
## 3. Wavelet transform (row-wise)
## -------------------------

# Perform wavelet decomposition for each row
W <- lapply(1:n, function(i)
  wd(Y[i, ])
)

# Extract all wavelet coefficients into a matrix
Wmat <-  do.call(rbind, lapply(1:length(W), function(i) W[[i]]$D))

## -------------------------
## 4. Correlation: wavelet domain
## -------------------------

cor_wavelet <- cor(Wmat)

wave_corr_mean <- mean(abs(cor_wavelet[upper.tri(cor_wavelet)]))

## -------------------------
## 5. Scale-wise correlations
## -------------------------

nlevels <- W[[1]]$nlevels

scale_corr <- numeric(nlevels)

for (j in 1:(S-1)) {
  Dj <- t(sapply(W, function(w)
    accessD(w, level = j)
  ))

  Cj <- cor(Dj)
  scale_corr[j] <- mean(abs(Cj[upper.tri(Cj)]))
}

## -------------------------
## 6. Results
## -------------------------

cat("Average |correlation|\n")
cat("----------------------\n")
cat("Original domain :", orig_corr_mean, "\n")
cat("Wavelet domain  :", wave_corr_mean, "\n\n")

cat("Scale-wise wavelet correlations:\n")
print(scale_corr)

plot_corr <- function(C, main) {
  C=t(C)
  image(
    1:nrow(C), 1:ncol(C), C[nrow(C):1, ],
    col = colorRampPalette(c("blue", "white", "red"))(200),
    zlim = c(-1, 1),
    axes = FALSE,
    main = main
  )
  box()
}


par(mfrow = c(2, 2))
plot_corr(cor_time,    "Original domain correlation")
plot_corr(cor_wavelet, "Wavelet domain correlation")



hist(c(cor_time), nclass = 200 ,main= "Original domain correlation")

hist(c(cor_wavelet), nclass = 200 ,main= "Wavelet domain correlation")
par(mfrow = c(1, 1))


pdf(paste0 (path, "/plot/correlation_gp_wavelet.png") , width = 11.69, height = 11.69)  
par(mfrow = c(2, 2))
plot_corr(cor_time,    "Original domain correlation")
plot_corr(cor_wavelet, "Wavelet domain correlation")



hist(c(cor_time), nclass = 200 ,main= "Original domain correlation")

hist(c(cor_wavelet), nclass = 200 ,main= "Wavelet domain correlation")
par(mfrow = c(1, 1))

dev.off()