library(fsusieR)
load("/home/wdenault/fsusi_simu/sim/genotypes.Rdata")
geno_info = readRDS("/home/wdenault/fsusi_simu/sim/Yuqi_data/geno_list_MWE.rds")
genotype <- t(geno_info$geno[1:2000,1:100])
library(gplots)

if(file.exists("/home/wdenault/fsusi_simu/correlated_sim/long_range_512_sd1.RData")){
  load("/home/wdenault/fsusi_simu/correlated_sim/long_range_512_sd1.RData")

}else{
  res <-list()
}


library(fracdiff)

S=9

H=0.75

target_pve <- 0.10
for (o  in (length(res)+1):10000) {
  set.seed(o)
  L <- sample(1:10, size=1)#actual number of effect
  lf <-  list()
  for(l in 1:L){
    lf[[l]] <- simu_IBSS_per_level(lev_res=7)$sim_func #functional effect for effect l
  }


  tt <- sample(0:4,1)
  G <- genotypes

  if( length(which(apply(G,2,var)==0))>0){
    G <- G[,-which(apply(G,2,var)==0)]
  }
  # G <- matrix( rnorm(100*300), nrow = 100)
  true_pos <- sample( 1:ncol(G), L)

  Y_signal <- matrix(0, nrow = 100, ncol = 2^S)

  for (i in 1:100) {
    for (l in 1:L) {
      Y_signal[i, ] <- Y_signal[i, ] + lf[[l]] * G[i, true_pos[l]]
    }
  }
  Y_noise <- matrix(0, nrow = 100, ncol = 2^S)

  for (i in 1:100) {
    eps <- fracdiff.sim(n = 2^S, d = H - 0.5)$series
    Y_noise[i, ] <- eps / sd(eps)
  }

  var_signal <- var(as.vector(Y_signal))
  var_noise  <- var(as.vector(Y_noise))

  signal_scale <- sqrt(
    (target_pve / (1 - target_pve)) * (var_noise / var_signal)
  )
  Y <- signal_scale * Y_signal + Y_noise



  m1 <- susiF(Y=Y, X=G,L=20 ,L_start=11 ,nullweight= 10,  prior="mixture_normal", cal_obj =FALSE,  maxit=10)
  m2 <- susiF(Y=Y, X=G,L=20,L_start=11 ,nullweight= 10,
              prior="mixture_normal_per_scale"  , maxit=10)

  cal_purity <- function(l_cs,X){
    tt <- list()
    for (k in 1:length(l_cs)){
      if(length(unlist(l_cs[[k]]))==1 ){
        tt[[k]] <- 1
      }else{
        x <-abs( cor(X[,unlist(l_cs[[k]]   ) ]))


        tt[[k]] <-  min( x[col(x) != row(x)])
      }
    }
    return( mean( unlist(tt )))
  }

  out <- c( length(m1$cs), #number of CS
            length(which(true_pos%in% do.call(c, m1$cs))), #number of effect found
            Reduce("+",sapply(1:length(m1$cs), function(k)
              ifelse( length(which(true_pos%in%m1$cs[[k]] ))==0, 1,0)
            )
            ),#number of CS without any effect
            cal_purity(m1$cs, X=as.matrix(G)),#mean purity
            mean(sapply( m1$cs, length)), #CS size
            length(m2$cs),
            length(which(true_pos%in% do.call(c, m2$cs))),
            Reduce("+",sapply(1:length(m2$cs), function(k)
              ifelse( length(which(true_pos%in%m2$cs[[k]] ))==0, 1,0)
            )
            ),#number of CS without any effect
            cal_purity(m2$cs, X=as.matrix(G)),#mean purity
            mean(sapply( m2$cs, length)), #CS size
            L,tt)
  res[[o]] <-unlist(out)


  save(res, file="/home/wdenault/fsusi_simu/correlated_sim/long_range_512_sd1.RData")
}

