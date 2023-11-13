library(susiF.alpha)
library(susieR)
load("/home/wdenault/fsusi_simu/sim/genotypes.Rdata")
geno_info = readRDS("/home/wdenault/fsusi_simu/sim/Yuqi_data/geno_list_MWE.rds")
data(N3finemapping)
X <- N3finemapping$X
N=50
genotype <-X[1:N,]

idx <- which( apply( genotype,2, var ) <1e-15)
genotype <- genotype [, -idx]
library(gplots)#

if(file.exists("/home/wdenault/fsusi_simu/sim3/overlap_check_gaussian.RData")){
  load("/home/wdenault/fsusi_simu/sim3/overlap_check_gaussian.RData")

}else{
  res <-list()
}

Rtrue <- cor (genotype )
for (o  in (length(res)+1):10000) {

  L <- sample(1, size =1)#actual number of effect
  lf <-  list()
  for(l in 1:L){
    lf[[l]] <- 0.5*simu_IBSS_per_level(lev_res=7)$sim_func #functional effect for effect l
  }


  genotype <-X[sample (1:nrow( X), size = N, replace =FALSE),]
  idx <- which( apply( genotype,2, var ) <1e-15)
  if(length( idx)>0 ){

    genotype <- genotype [, -idx]

  }


  tt <- sample(0:4,1)
  G <- genotype

  if( length(which(apply(G,2,var)==0))>0){
    G <- G[,-which(apply(G,2,var)==0)]
  }
  # G <- matrix( rnorm(N*300), nrow = N)
  true_pos <- sample( 1:ncol(G), L)

  Y <- matrix(0 , ncol=  2^7 , nrow = N)
  for ( i in 1:N){
    for ( l in 1:L){
      Y[i,] <- Y[i,]+ lf[[l]]*G[i,true_pos[[l]]]
    }
  }


  Y <-Y+matrix(rnorm((2^7)*N ,sd=sd(c(Y))/sqrt(0.5)), nrow = N)

  PCA <- svd(Y)
  m1 <-susiF(Y=Y, X=G,L=2  ,nullweight=10,  maxit=10)
  ### buidling zscore


  m2 <-susie(X=G,
             y=PCA$u[,1],
             L=5
  )



  lf <-  list()
  for(l in 1:L){
    lf[[l]] <- 0.5*simu_IBSS_per_level(lev_res=7)$sim_func #functional effect for effect l
  }


  Y2 <- matrix(0 , ncol=  2^7 , nrow = N)
  for ( i in 1:N){
    for ( l in 1:L){
      Y2[i,] <- Y2[i,]+ lf[[l]]*G[i,true_pos[[l]]]
    }
  }


  Y2 <-Y2+matrix(rnorm((2^7)*N ,sd=sd(c(Y))/sqrt(0.5)), nrow = N)

  PCA <- svd(Y2)
  m12 <-susiF(Y=Y2, X=G,L=2  ,nullweight=10,  maxit=10)
  ### buidling zscore


  m22 <-susie(X=G,
             y=PCA$u[,1],
             L=2
  )


  out <-  list(   is_overlap_susif  =       ifelse(length( which(m1$cs[[1]] %in% m12$cs[[1]] )) >0,
                                               1,
                                               0),
                  m1$pip,
                  susie_rss_pip= m2$pip,
                  is_overlap_susie=  ifelse(length( which(m2$sets$cs$L1 %in% m22$sets$cs$L1 )) >0,
                                            1,
                                            0) ,
                  susiF_pip = m1$pip, susie_rss_pip= m2$pip,
                  susiF_cs= m1$cs,
                  susie_cs= m2$sets,
                  fsusie_purity = cal_purity(m1$cs,G  ),
                  true_pos=true_pos)

  res[[o]] <- out
  save(res, file="/home/wdenault/fsusi_simu/sim3/overlap_check_gaussian.RData")

}
