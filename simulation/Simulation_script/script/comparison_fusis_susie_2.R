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

if(file.exists("/home/wdenault/fsusi_simu/sim3/comparison_susie_fusie_distdecay_sd1.RData")){
  load("/home/wdenault/fsusi_simu/sim3/comparison_susie_fusie_distdecay_sd1.RData")

}else{
  res <-list()
}

Rtrue <- cor (genotype )
for (o  in (length(res)+1):10000) {
  L <- sample(1:5, size =1)#actual number of effect
  lf <-  list()
  R=2^7
  for(l in 1:L){

    block_cov_generate = function(corr=0.7, v=1){
      R = dim(corr)[1]
      cov = diag(sqrt(v),R) %*% corr %*% diag(sqrt(v),R)
      return(cov)
    }

    decay_corr_generate = function(R =2^7,# # sites
                                   num_block = 10,# #islands within the region
                                   off_diag_corr = 0.7,
                                   decay =0.9#
    ){
      corr_b = matrix(0,nrow=R,ncol=R)
      region_mat = matrix(off_diag_corr,R/num_block,R/num_block)
      for ( i in 1: nrow(region_mat)){
        for ( j in  1:ncol(region_mat)){
          tt <- region_mat[i,j]
          region_mat[i,j] <- tt*decay ^( abs(i-j))
          region_mat[j,i] <- tt*decay ^( abs(i-j))

        }
      }
      for(i in 1:num_block){
        corr_b[seq(((i-1)*(R/num_block)+1),i*(R/num_block)), seq(((i-1)*(R/num_block)+1),i*(R/num_block))] = region_mat
      }

      diag(corr_b) = 1
      return(corr_b)
    }
    corr <-    decay_corr_generate()
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


  Y <-Y+matrix(rnorm((2^7)*N ,sd=sd(c(Y))/sqrt(1)), nrow = N)


  sigmoid <- function( x)
  {
    out <- 1/(1+exp(-x))
    return(out)
  }
  Y <- sigmoid(Y)
  PCA <- svd(Y)
  m1 <-susiF(Y=Y, X=G,L=5  ,nullweight=10,  maxit=10)
  ### buidling zscore


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



  out <-  list( susiF_pip = m1$pip, susie_rss_pip= m2$pip,
                susiF_cs= m1$cs,
                susie_cs= m2$sets,
                fsusie_purity = cal_purity(m1$cs,X  ),
                true_pos=true_pos)

  res[[o]] <- out
  save(res, file="/home/wdenault/fsusi_simu/sim3/comparison_susie_fusie_distdecay_sd1.RData")
  print(o)
}

