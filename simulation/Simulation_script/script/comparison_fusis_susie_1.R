library(susiF.alpha)
library(susieR)
load("/home/wdenault/fsusi_simu/sim/genotypes.Rdata")
geno_info = readRDS("/home/wdenault/fsusi_simu/sim/Yuqi_data/geno_list_MWE.rds")
data(N3finemapping)
X <- N3finemapping$X
N=50
genotype <-X[1:N,1:300]

idx <- which( apply( genotype,2, var ) <1e-15)
if ( length( idx )>0){

  genotype <- genotype [, -idx]
}
library(gplots)#

if(file.exists("/home/wdenault/fsusi_simu/sim3/comparison_susie_fusie_128_sd1.RData")){
  load("/home/wdenault/fsusi_simu/sim3/comparison_susie_fusie_128_sd1.RData")

}else{
  res <-list()
}

Rtrue <- cor (genotype )
for (o  in (length(res)+1):10000) {
  set.seed (o)
  L <- sample(1:5, size =1)#actual number of effect
  lf <-  list()
  for(l in 1:L){
    lf[[l]] <- simu_IBSS_per_level(lev_res=7)$sim_func #functional effect for effect l
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


  Y <-Y+matrix(rnorm((2^7)*N ,sd=sd(c(Y))/sqrt(1)), nrow = N)

  PCA <- svd(Y)
  m1 <-susiF(Y=Y, X=G,L=5  ,nullweight=10,  maxit=10)
  ### buidling zscore
  m11 <-susiF(Y=Y, X=G,L=5  ,nullweight=10,  maxit=10,
              prior="mixture_normal",
              post_processing="none")

  m2 <-susie(X=G,
             y=PCA$u[,1],
             L=5
  )
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




  out <-  list( susiF_pip = m1$pip,
                susiF_sp_pip = m11$pip,

                susie_rss_pip= m2$pip,
                susiF_cs= m1$cs,
                susiF_sp_cs = m11$cs,
                susie_cs= m2$sets,

                fsusie_purity = cal_purity(m1$cs,G ),
                fsusie_sp_purity = cal_purity(m11$cs,G ),
                true_pos=true_pos)

  res[[o]] <- out
  save(res, file="/home/wdenault/fsusi_simu/sim3/comparison_susie_fusie_128_sd1.RData")

}
