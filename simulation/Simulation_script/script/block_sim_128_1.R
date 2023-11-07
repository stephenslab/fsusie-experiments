library(susiF.alpha)
load("/home/wdenault/fsusi_simu/sim/genotypes.Rdata")
geno_info = readRDS("/home/wdenault/fsusi_simu/sim/Yuqi_data/geno_list_MWE.rds")
genotype <- t(geno_info$geno[1:2000,1:100])
library(gplots)
if(file.exists("/home/wdenault/fsusi_simu/sim3/block_L_accuracy_128_sd1.RData")){
  load("/home/wdenault/fsusi_simu/sim3/block_L_accuracy_128_sd1.RData")

}else{
  res <-list()
}
for (o  in (length(res)+1):10000) {

  L <- sample(1:20, size=1)#actual number of effect
  lf <-  list()
  R=2^7
  for(l in 1:L){

    block_cov_generate = function(corr=0.8, v=1){
      R = dim(corr)[1]
      cov = diag(sqrt(v),R) %*% corr %*% diag(sqrt(v),R)
      return(cov)
    }

    block_corr_generate = function(R =2^7,# # sites
                                   num_block = 1,# #islands within the region
                                   off_diag_corr = 0.8 # correlation within the island
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
  G <- genotypes

  if( length(which(apply(G,2,var)==0))>0){
    G <- G[,-which(apply(G,2,var)==0)]
  }
  # G <- matrix( rnorm(100*300), nrow = 100)
  true_pos <- sample( 1:ncol(G), L)

  Y <- matrix(0 , ncol=  2^7 , nrow = 100)
  for ( i in 1:100){
    for ( l in 1:L){
      Y[i,] <- Y[i,]+ lf[[l]]*G[i,true_pos[[l]]]
    }
  }


  Y <-Y+matrix(rnorm((2^7)*100 ,sd=sd(c(Y))), nrow = 100)


  sigmoid <- function( x)
  {
    out <- 1/(1+exp(-x))
    return(out)
  }
  Y <- sigmoid(Y)
  m1 <- susiF(Y=Y, X=G,L=20 ,L_start=11 ,nullweight=10,  prior="mixture_normal", cal_obj =FALSE,  maxit=10)
  m2 <- susiF(Y=Y, X=G,L=20,L_start=11 ,nullweight=10 , maxit=10)
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


  save(res, file="/home/wdenault/fsusi_simu/sim3/block_L_accuracy_128_sd1.RData")
}

