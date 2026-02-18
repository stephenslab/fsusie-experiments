library(mvsusieR) 
genotypes = N3finemapping$X[ sample(1:nrow(N3finemapping$X), size=100),]
set.seed(1)

s=7

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
prior <- create_mixture_prior(R = ncol(  Y ))
m1= mvsusie(X=G, Y= Y , prior_variance = prior)
rt_mvsusie= proc.time()-pt
 
 