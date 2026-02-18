# check_lfsr_calibration_functional.R
library(fsusieR)
rm(list = ls())
set.seed(1)
source("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/additional_work_for_revision/script/RMSE_comparison.R", echo = F)
source("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/additional_work_for_revision/script/RMSE_comparison_block.R", echo = F)

path <- getwd()
#Simulation to hcek on Bump functoion calibartion of hmm lfsr
#two modality SNP affect only the first on



# -------------------------
# Parameters (tune as needed)
# -------------------------
N <- 100            # samples
P <-1           # SNPs
Tlen <- 128       # functional domain length (2^6)
nsim <-50       # number of replicates (300 is a good compromise)
L_fit <- 1


# bump specification (location and amplitude)
bump_center <- round(Tlen * 0.5)
bump_width <- round(Tlen * 0.12)    # half-width of bump
bump_amp <- 1.5                     # amplitude of bump

# utility: build bump vector length Tlen
build_bump <- function(Tlen, center, width, amp) {
  x <- seq_len(Tlen)
  w <- width
  bump <- amp * exp(-((x - center)^2) / (2 * (w/2)^2))
  return(bump)
}
true_bump <- build_bump(Tlen, bump_center, bump_width, bump_amp)
true_bump[which(true_bump<1e-2)]=0
# Option: you can set some domain points to exactly zero to test null calibration outside bump
is_nonzero_point <- which(true_bump != 0)

# storage
# For pointwise calibration we'll collect, across replicates, the lfsr vector when
# there is a causal SNP. We treat points outside the bump as "null" (true effect ~ 0).
pointwise_calls <- matrix(0L, nrow=nsim, ncol=Tlen)  # 1 if lfsr < thresh (we'll evaluate per threshold)
trait_min_lfsr <- numeric(nsim)  # min lfsr across domain per replicate
cs_contains_causal <- logical(nsim)

result=list()
o=1
for ( sds in c(0.1,.5,1,2,4)){
  
  
  # Optionally allow null-only replicates (no causal SNP) to assess behavior under pure null:
  # Here we simulate the causal SNP present always, but the bump amplitude can be set to 0 to test null.
  for (rep in 1:nsim) {
    # 1) Simulate genotype
    
    
    set.seed(rep)
    off_diag_corr = 0.99
    decay =0.9
    
    region_mat= matrix(off_diag_corr, ncol=  128, nrow=  128)
    for ( i in 1: nrow(region_mat)){
      for ( j in  1:ncol(region_mat)){
        tt <- region_mat[i,j]
        region_mat[i,j] <- tt*decay ^( abs(i-j))
        region_mat[j,i] <- tt*decay ^( abs(i-j))
        
      }
    }
    
    effec =dae::rmvnorm(mean = rep(0, 128),V = region_mat) 
    # bump specification (location and amplitude)
    bump_center <- round(Tlen * 0.5)
    bump_width <- round(Tlen * 0.12)    # half-width of bump
    bump_amp <- 1                    # amplitude of bump
    
    # utility: build bump vector length Tlen
    build_bump <- function(Tlen, center, width, amp) {
      x <- seq_len(Tlen)
      w <- width
      bump <- amp * exp(-((x - center)^2) / (2 * (w/2)^2))
      return(bump)
    }
    true_bump <- build_bump(Tlen, bump_center, bump_width, bump_amp)
    true_bump[which(true_bump<1e-2)]=0
    true_bump=  true_bump*effec
    # 1) Simulate genotype
    G <- matrix(sample(0:2, size=N*P, replace=TRUE), nrow=N, ncol=P)
    causal_snp <- sample(1:P, 1)
    
    
    # 2) Simulate functional phenotype Y_f (N x Tlen): Gaussian noise + genotype*bump at causal SNP
    Yf <- matrix(rnorm(N * Tlen, sd=sds  ), nrow=N, ncol=Tlen)
    Yf2 <- matrix(rnorm(N * Tlen, sd=sds ), nrow=N, ncol=Tlen)
    # add genotype-mediated bump
    for (i in 1:N) {
      Yf[i, ] <- Yf[i, ] + G[i, causal_snp] * true_bump
    }
    Y <- list(Y_f = list(Yf, Yf2), Y_u = NULL)   # single functional trait
    
    # pos argument: domain positions (needed by multfsusie)
    
    
    # 3) Fit multfsusie
    Y = Yf
    X = G
    mfit <-   susiF(Y = Yf, X = G , L = L_fit, post_processing ="smash", verbose = FALSE)
    mfit2 <-   susiF(Y = Yf, X = G , L = L_fit, post_processing ="TI", verbose = FALSE)
    
    
    est_naive=unlist(lapply(1: ncol( Yf),  function (i){
      
      
      summary(lm(Yf[,i] ~ X))$coefficient[2,1]
    }))
    
    
    # 4) Extract pointwise functional LFSR.
    #    mfit$lfsr is a list of effects; each mfit$lfsr[[e]]$est_lfsr_functional[[1]] is a vector length Tlen
    # Aggregate across effects by taking minimum (conservative): for each domain point take min LFSR across effects
    
    result[[o]]= list( cb1=mfit$cred_band[[1]]  ,
                       f1=mfit$fitted_func[[1]]  ,
                       f2=mfit2$fitted_func[[1]]  ,
                       est_naive=  est_naive,
                       sd=sds ,
                       true_bump=true_bump )
    print(rep)
    o=o+1
  }
}
rmse=function(x,y){
  sqrt(mean((x-y)^2))
}


est_IS=do.call( rbind, lapply( 1: length(result), function(i)  rmse (result[[i]]$f2,result[[i]]$true_bump)))

est_smash=do.call( rbind, lapply( 1: length(result), function(i)  rmse (result[[i]]$f1,result[[i]]$true_bump)))

est_ols=do.call( rbind, lapply( 1: length(result), function(i)  rmse (result[[i]]$ est_naive,result[[i]]$true_bump)))


sds= do.call( rbind, lapply( 1: length(result), function(i)  result[[i]]$ sd ))

library(ggplot2)
df= data.frame(sds=factor(rep(sds,3)),
               rmse= c(est_ols,est_smash,est_IS),
               method =factor(rep(c("SNP-CpG association tests", "fSuSiE SPSP", "fSuSiE ISP"), each=length(sds))))

colors= c(  "#D41159","#FFC20A" ,"darkblue")#c(  "#FFC20A","#D41159" )
P3= ggplot(df, aes(x=sds, y=rmse, col=method))+
  geom_boxplot() +
  xlab("Noise level")+
  ylab("RMSE")+
  scale_color_manual(values=colors)+
  ggtitle("WGBS decay")+
  theme_linedraw()
P3
pdf(paste0 (path, "/plot/RMSE_comp_decay.pdf") , width = 10.27, height = 8.27)  
P3

dev.off()


library(gridExtra)

pdf(paste0 (path, "/plot/RMSE_sim.pdf") , width = 10.27, height = 5.27)  

library(patchwork)
 P2= P2+ggtitle("wavelet simulations"
                )
 P3= P3+ggtitle("WGBS block"
 )
 P1= P1+ggtitle("WGBS decay"
 )
(  P2 + P3+P1) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom") &
  guides(color = guide_legend(title = NULL))

dev.off()
