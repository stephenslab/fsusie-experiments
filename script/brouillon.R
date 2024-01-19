effect <-   0.2*cos( (1:128)/128 * 3*pi )
effect[which(effect<0)]<- 0
effect[1:40]<- 0

plot( effect)
library(susiF.alpha)
library(susieR)
library(wavethresh)
set.seed(2)
data(N3finemapping)
X <- N3finemapping$X[1:120,]

set.seed(1)

obs <- list()

true_pos <- 350

table (X[,true_pos])
X[,true_pos] <- (X[,true_pos] - min(X[,true_pos] ))
table (X[,true_pos])
for (i in 1:120 ){
  obs[[i]] <- X[i,true_pos]*effect+ runif(length(effect),max=.5)
}

y <- do.call(rbind, obs)
y <- y[, 50:120]

susie_est <- list()

for (i in 1:ncol(y)){
  
  
  susie_est [[i]]<- susie( X=X,y=y[,i], L=1)$sets
}

susie_est
#48
#45
#44
#43
#42
#41
library(ggplot2)
#88
set.seed(12345)
pos <- cumsum( runif(n=length(50:120)))

idx <-  sample (1:length(50:120), size=30)
idx <- idx[order(idx)]

plot(pos[idx],effect[50:120][idx])

susie_est <- list()
h=1
y <- do.call(rbind, obs)
y <- y[, 50:120]
for (i in idx){
  
  
  susie_est [[h]]<- susie( X=X,y=y[,i], L=1)$sets
  h=h+1
}

susie_est

colnames(X) <- 1:ncol(X)

tt <- susiF(y[,idx] , X, L=1 ,pos=pos[idx] )
tt$cs
plot(pos[idx],effect[50:120][idx])

lines(x=tt$outing_grid , y=tt$fitted_func[[1]])

lines(x=tt$outing_grid , y=tt$cred_band[[1]][1,])
lines(x=tt$outing_grid , y=tt$cred_band[[1]][2,])

p <- ggplot()
to_remove <- which(X[,true_pos]==0.39892578125)
tl <- list()
h <- 1
for ( i in idx){
  
  tl[[h]]  <- data.frame(y =  y[-to_remove,i],
                         x =rep( pos[i ] , nrow(y)-1),
                         Genotype=  X[-to_remove ,true_pos]
  )
  h <- h+1
  
}
df <- do.call(rbind, tl) 


table(df$Genotype)
lst_pos <-list()
for ( i in 1:30){
  lst_pos[[i]] <- 1*rep( unique(df$x)[i],3)+ c(0,0.05,0.1)
}

 
myplot <- boxplot(y ~ as.factor(Genotype)*x , data=df  , 
                  boxwex=0.2 ,
                  ylab="", 
                  yaxt="n",
                  xlab= "",
                  at= do.call( c,lst_pos),
                  col=c("tomato","slateblue1" , "green2") ,  
                  xaxt="n")
lines(x= .5+1:length(50:120)*(35/(120-50)), y=effect[50:120]+0.25,
      col="slateblue1", lwd=2)
lines(x=.5+1:length(50:120)*(35/(120-50)), y=2*effect[50:120]+0.25,
      col="green2", lwd=2)
lines(x=.5+1:length(50:120)*(35/(120-50)), y=0*effect[50:120]+0.25,
      col="tomato", lwd=2)




library(ggplot2)
ggplot(df, aes(x=x, y=y, col=as.factor(Genotype)))+
  geom_point()
