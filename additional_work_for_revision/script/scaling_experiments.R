library(fsusieR)
load("/home/wdenault/fsusi_simu/sim/genotypes.Rdata")
geno_info = readRDS("/home/wdenault/fsusi_simu/sim/Yuqi_data/geno_list_MWE.rds")
genotype <- t(geno_info$geno[1:2000,1:100])
library(gplots)

library(mvsusieR)
library(susieR)
if(file.exists("/home/wdenault/fsusi_simu/correlated_sim/run_time.RData")){
  load("/home/wdenault/fsusi_simu/correlated_sim/run_time.RData")

}else{
  res <-list()
}

o=1

for (s  in 5:10) {



for ( k in 1:10){

  L <- 3 #actual number of effect
  lf <-  list()
  for(l in 1:L){
    lf[[l]] <- simu_IBSS_per_level(lev_res=s)$sim_func #functional effect for effect l
  }


  tt <- sample(0:4,1)
  G <- genotypes

  if( length(which(apply(G,2,var)==0))>0){
    G <- G[,-which(apply(G,2,var)==0)]
  }
  # G <- matrix( rnorm(100*300), nrow = 100)
  true_pos <- sample( 1:ncol(G), L)

  Y <- matrix(0 , ncol=  2^s , nrow = 100)
  for ( i in 1:100){
    for ( l in 1:L){
      Y[i,] <- Y[i,]+ lf[[l]]*G[i,true_pos[[l]]]
    }


    Y[i,] <- Y[i,]
  }

  Y <- Y   +matrix(rnorm((2^s)*100 ,sd=sd(c(Y))), nrow = 100)


  pt = proc.time()
  m1 <- susiF(Y=Y, X=G,L=20 ,L_start=5 ,nullweight=10,tol = 1e-6,
              prior="mixture_normal", cal_obj =FALSE,  maxit=10, post_processing = "none")
  rt_fsusie= proc.time()-pt

  pt_inner=0*rt_fsusie
  pt = proc.time()
  for (j in 1:ncol(Y) ) {
    pt0 = proc.time()
    m1 <- susie(y=Y[,j], X=G)
    ptemp=    proc.time()-pt0
    pt_inner=pt_inner+ ptemp
  }
  rt_susie= proc.time()-pt

  #pt = proc.time()
  #prior <- create_mixture_prior(R = ncol(  Y ))
  #m1= mvsusie(X=G, Y= Y , prior_variance = prior)
  #rt_mvsusie= proc.time()-pt

  out <-  list(rt_fsusie=rt_fsusie,
               rt_susie= rt_susie,
               pt_inner= pt_inner,
               s=s
               #rt_mvsusie=rt_mvsusie
               )
  res[[o]] <- (out)


  save(res, file="/home/wdenault/fsusi_simu/correlated_sim/run_time.RData")
  o=o+1
  }
}


