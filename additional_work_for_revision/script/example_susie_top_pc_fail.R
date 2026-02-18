 
library(fsusieR)
library(reshape2)
library(ggplot2)
library(cowplot)
 library(susieR)
set.seed(1)
  
n <- 100
m <- 32
p <- 12

cs_colors <- c("dodgerblue","darkorange","red")
 
maf <- 0.05 + 0.45*runif(p)
 
snpids <- paste0("SNP-",1:p)
cpgids <- paste0("CpG-",1:m)
 
X <- (runif(n*p) < maf) +
  (runif(n*p) < maf)
X <- matrix(X,n,p,byrow = TRUE)
storage.mode(X) <- "double"
X[,4] <- X[,3] + 0.03*rnorm(n)
colnames(X) <- snpids
  
F <- matrix(0,p,m)
F[1,9:16] <- 2.3
F[9,9:16] <- (-2.3)
F[3,25:32] <- 2
rownames(F) <- snpids
colnames(F) <- cpgids
 
E <- matrix(3*rnorm(n*m),n,m)
Y <- X %*% F
Y <- Y + E
baseline <- min(Y)
Y <- Y - baseline
 
pdat <- melt(Y)
x    <- X[,3]
pdat <- data.frame(cpg        = rep(1:m,each = n) +
                     runif(m*n,min = -0.2,max = 0.2),
                   meth_level = as.vector(Y),
                   geno       = factor(x))
pdat_lines <- data.frame(cpg = rep(1:m,times = 3),
                         geno = factor(rep(0:2,each = m)),
                         meth_level = rep(c(rep(0,8),rep(-0.75,8),rep(0,16)),
                                          times = 3))
pdat_lines$meth_level <- pdat_lines$meth_level - baseline
rows <- which(with(pdat_lines,geno == 1 & cpg > 24))
pdat_lines[rows,"meth_level"] <- pdat_lines[rows,"meth_level"] + 2
rows <- which(with(pdat_lines,geno == 2 & cpg > 24))
pdat_lines[rows,"meth_level"] <- pdat_lines[rows,"meth_level"] + 4
p1 <- ggplot(pdat, aes(x = cpg, y = meth_level, color = geno)) +
  geom_point(shape = 20, size = 1.25) +
  scale_x_continuous(breaks = c(0, seq(4, 32, 4))) +
  scale_color_manual(values = c("darkblue", "limegreen", "darkorange")) +
  geom_line(data = pdat_lines, aes(x = cpg, y = meth_level, group = geno), 
            color = "black", size = 1.9) +  # Black line for contrast
  geom_line(data = pdat_lines, size = 1.25) +  # Original colored lines
  labs(x = "CpG", y = "methylation level", color = "genotype") +
  theme_cowplot(font_size = 11)
print(p1)
  
 
fit_TI <- susiF(Y,X,L = 3,filter_cs = FALSE,prior = "mixture_normal" )
  
fit_TI$cs
 
out <- lapply(1:3,function (i) get_fitted_effect(fit_TI,l = i,cred_band = TRUE,
			                                     alpha = 0.05))
effect_plot <- function (i) {
  pdat <- data.frame(cpg      = 1:m,
                     estimate = out[[i]]$effect,
					 lower    = out[[i]]$cred_band["low",],
					 upper    = out[[i]]$cred_band["up",])
  rows <- with(pdat,which(lower > 0 | upper < 0))
  pdat2 <- data.frame(cpg = rows,estimate = 0)
  return(ggplot(pdat,aes(x = cpg,y = estimate,ymin = lower,ymax = upper)) +
         geom_point(color = cs_colors[i],size = 1) +
         geom_errorbar(color = cs_colors[i],linewidth = 0.4) +
         geom_hline(yintercept = 0,color = "black",linewidth = 0.4) +
		 geom_point(data = pdat2,mapping = aes(x = cpg,y = estimate),
		            shape = 20,color = "black",size = 1.5,
					inherit.aes = FALSE) +
         scale_x_continuous(breaks = c(0,seq(4,32,4))) +
         labs(x = "CpG",y = "change",title = paste0("CS",i)) +
		 theme_cowplot(font_size = 11))
}
plot_grid(effect_plot(1),
          effect_plot(2),
   	      effect_plot(3),
		  nrow = 3,ncol = 1) 
 




pca_res= prcomp(Y, center = TRUE, scale. = TRUE)



res_susie_PC1=susie(y=pca_res$x[,1],X)
res_susie_PC1$sets

res_susie_PC2=susie(y=pca_res$x[,2],X)
res_susie_PC2$sets
 