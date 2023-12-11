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
library(gplots)
if(file.exists("/home/wdenault/fsusi_simu/sim3/overlap_check_block_sd1.RData")){
  load("/home/wdenault/fsusi_simu/sim3/overlap_check_block_sd1.RData")

}else{
  res <-list()
}
Rtrue <- cor (genotype )
for (o  in (length(res)+1):10000) {

  L <- sample(1:1, size=1)#actual number of effect
  lf <-  list()


  genotype <-X[sample (1:nrow( X), size = N, replace =FALSE),]
  idx <- which( apply( genotype,2, var ) <1e-15)
  if(length( idx)>0 ){

    genotype <- genotype [, -idx]

  }
  R=2^7
  for(l in 1:L){

    block_cov_generate = function(corr=0.8, v=1){
      R = dim(corr)[1]
      cov = diag(sqrt(v),R) %*% corr %*% diag(sqrt(v),R)
      return(cov)
    }

    block_corr_generate = function(R =2^7,# # sites
                                   num_block = 10,# #islands within the region
                                   off_diag_corr = 0.5 # correlation within the island
                                   # v = 0.001, # residual variance(before logit transformation)
                                   # pve = 0.01, # heritability par
                                   # pi_0_type = "beta"
    ){
      corr_b = matrix(0,nrow=R,ncol=R)
      region_mat = matrix(off_diag_corr,R/num_block,R/num_block)
      for(i in 1:num_block){
        corr_b[seq(((i-1)*(R/num_block)+1),i*(R/num_block)), seq(((i-1)*(R/num_block)+1),i*(R/num_block))] = region_mat
      }

      diag(corr_b) = 1
      return(corr_b)
    }
    corr <-  block_corr_generate()
    cov = block_cov_generate(corr =   corr, v = 0.1)


    lf[[l]] <- dae::rmvnorm(mean = rep(0,R),V = cov) #functional effect for effect l
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


  Y <-Y+matrix(rnorm((2^7)*N ,sd=sd(c(Y))), nrow = N)


  sigmoid <- function( x)
  {
    out <- 1/(1+exp(-x))
    return(out)
  }
  Y <- sigmoid(Y)

  PCA <- svd(Y)
  m1 <-susiF(Y=Y, X=G,L=5  ,nullweight=10,  maxit=10)
  ### buidling zscore

  m01 <-susiF(Y=Y, X=G,L=2  ,nullweight=10,  maxit=10,
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





  for(l in 1:L){

    block_cov_generate = function(corr=0.8, v=1){
      R = dim(corr)[1]
      cov = diag(sqrt(v),R) %*% corr %*% diag(sqrt(v),R)
      return(cov)
    }

    block_corr_generate = function(R =2^7,# # sites
                                   num_block = 10,# #islands within the region
                                   off_diag_corr = 0.5 # correlation within the island
                                   # v = 0.001, # residual variance(before logit transformation)
                                   # pve = 0.01, # heritability par
                                   # pi_0_type = "beta"
    ){
      corr_b = matrix(0,nrow=R,ncol=R)
      region_mat = matrix(off_diag_corr,R/num_block,R/num_block)
      for(i in 1:num_block){
        corr_b[seq(((i-1)*(R/num_block)+1),i*(R/num_block)), seq(((i-1)*(R/num_block)+1),i*(R/num_block))] = region_mat
      }

      diag(corr_b) = 1
      return(corr_b)
    }
    corr <-  block_corr_generate()
    cov = block_cov_generate(corr =   corr, v = 0.1)


    lf[[l]] <- dae::rmvnorm(mean = rep(0,R),V = cov) #functional effect for effect l
  }



  Y2 <- matrix(0 , ncol=  2^7 , nrow = N)
  for ( i in 1:N){
    for ( l in 1:L){
      Y2[i,] <- Y2[i,]+ lf[[l]]*G[i,true_pos[[l]]]
    }
  }


  Y2 <-Y2+matrix(rnorm((2^7)*N ,sd=sd(c(Y))/sqrt(1)), nrow = N)


  sigmoid <- function( x)
  {
    out <- 1/(1+exp(-x))
    return(out)
  }
  Y2 <- sigmoid(Y2)
  PCA <- svd(Y2)
  m12 <-susiF(Y=Y2, X=G,L=2  ,nullweight=10,  maxit=10,
              post_processing="none")
  ### buidling zscore
  m11 <-susiF(Y=Y, X=G,L=2  ,nullweight=10,  maxit=10,
              prior="mixture_normal",
              post_processing="none")

  ### buidling zscore


  m22 <-susie(X=G,
              y=PCA$u[,1],
              L=2
  )


  out <-  list(   is_overlap_susif  =       ifelse(length( which(m1$cs[[1]] %in% m12$cs[[1]] )) >0,
                                                   1,
                                                   0),
                  m1$pip,

                  is_overlap_susif_sp  =       ifelse(length( which(m11$cs[[1]] %in% m01$cs[[1]] )) >0,
                                                      1,
                                                      0),

                  susie_rss_pip= m2$pip,
                  is_overlap_susie=  ifelse(length( which(m2$sets$cs$L1 %in% m22$sets$cs$L1 )) >0,
                                            1,
                                            0) ,
                  susiF_pip = m1$pip,

                  susiF_sp_pip = m11$pip,

                  susie_rss_pip= m2$pip,
                  susiF_cs= m1$cs,
                  susiF_sp_cs= m11$cs,

                  susie_cs= m2$sets,
                  fsusie_purity = cal_purity(m1$cs,G  ),
                  true_pos=true_pos)
  res[[o]] <- out

  save(res, file="/home/wdenault/fsusi_simu/sim3/overlap_check_block_sd1.RData")
}

