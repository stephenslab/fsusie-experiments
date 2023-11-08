library(susiF.alpha)
library(susieR)
load("/home/wdenault/fsusi_simu/sim/genotypes.Rdata")


effect0 <-   1.2*cos( (1:128)/128 * 3*pi )


effect1 <- effect0
effect2 <- effect0


effect1[which(effect1<0)]<- 0
effect1[1:40]<- 0


effect2[which(effect2>0)]<- 0
effect2[80:length(effect0)]<- 0

library(susiF.alpha)
library(susieR)
data(N3finemapping)
X <- N3finemapping$X
genotype <-X[1:400,]


X <-genotype
cor_G<- cor(genotype)

library(gplots)#

if(file.exists("/home/wdenault/fsusi_simu/sim3/split_perf_400.RData")){
  load("/home/wdenault/fsusi_simu/sim3/split_perf_400.RData")

}else{
  res <-list()
}

for ( o in 1:10000){



  notOK <- TRUE
  while( notOK){
    pos_effect1 <- sample(1:ncol(X), size=1)
    if( length(which(cor_G[pos_effect1,]>0.5))>2){

      tt <- which(cor_G[pos_effect1,]>0.5)
      pos_effect2 <- sample(tt [-which(tt==pos_effect1)], size=1)

      notOK <-FALSE
    }
  }
  pos_effect1
  pos_effect2
  cor(genotype [, c(pos_effect1, pos_effect2)])[2,1]


  obs <- list()

  for (i in 1:nrow(X) ){
    obs[[i]] <- X[i,pos_effect1]*effect1+ X[i,pos_effect2]*effect2+ rnorm(length(effect0))
  }

  y <- do.call(rbind, obs)
  m1 <- susiF(Y=y, X=genotype,L=2)



  cor_SNP_effect <- cor(genotype [, c(pos_effect1, pos_effect2)])[2,1]
  cs_effect1 <- c()
  if(pos_effect1%in% m1$cs[[1]]){
    cs_effect1 <-c(cs_effect1 ,1)
  }
  if(length(m1$cs)>1){
    if(pos_effect1%in% m1$cs[[2]] ){
      cs_effect1 <-c(cs_effect1 ,2)
    }
  }


  cs_effect2 <- c()
  if(pos_effect1%in% m1$cs[[1]]){
    cs_effect2 <-c(cs_effect2 ,1)
  }


  if(length(m1$cs)>1){
    if(pos_effect1%in% m1$cs[[2]] ){
      cs_effect2 <-c(cs_effect2 ,2)
    }
  }
  out <- list( cs_effect1= cs_effect1,
               cs_effect2= cs_effect2,
               cor_SNP_effect = cor_SNP_effect,
               n_cs = length(m1$cs)
  )

  res[[o]] <- out
  save(res, file="/home/wdenault/fsusi_simu/sim3/split_perf_400.RData")
}
