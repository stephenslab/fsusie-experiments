rm(list=ls())
library(ashr)
library(wavethresh)
library(susiF.alpha)
library(ggplot2)
library(ggrepel)
set.seed(1)
#Example using curves simulated under the Mixture normal per scale prior
rsnr <- 0.2 #expected root signal noise ratio
N <- 100    #Number of individuals
P <- 10     #Number of covariates/SNP
pos1 <- 1   #Position of the causal covariate for effect 1
pos2 <- 5   #Position of the causal covariate for effect 2
lev_res <- 7#length of the molecular phenotype (2^lev_res)
f1 <-  simu_IBSS_per_level(lev_res )$sim_func#first effect
f2 <- simu_IBSS_per_level(lev_res )$sim_func #second effect
f2 <- simu_IBSS_per_level(lev_res )$sim_func #second effect


#### Function 1------
idx_lst <- gen_wavelet_indx(log2(128))


f1_list <- list()
wd1_list <- list()

xwd1 <-list()
ywd1 <-list()
for (s in 1:log2(length(f1))){
  tt <- rep( 0, 128)
  
  
  tt_wd <- wd(tt)
  tt_wd$D[ idx_lst[[s]]]  <- wd(f1)$D[ idx_lst[[s]]]
  
  wd1_list[ idx_lst[[s]]]  <- wd(f1)$D[ idx_lst[[s]]]
  xwd1[[s]] <- (idx_lst[[s]]-mean(idx_lst[[s]]))/length(idx_lst[[s]])
  ywd1[[s]] <- rep( s,length(idx_lst[[s]]))
  f1_list[[s]] <- wr(tt_wd)
}


df01 <- data.frame(func = c( f1 ),
                   x=  ( 1:128 )
                   
)


P01 <- ggplot( df01, aes(x=x, y=func))+
  geom_line(size=2,col="palegreen2")+
  theme_bw()+
  ylab("")+
  ggtitle("Effect to be estimated")+
  theme( )

df1 <- data.frame(func = c(  do.call(c, f1_list)),
                  x= rep( 1:128,  length(f1_list) ),
                  scale =  factor(
                    rep( 1:log2(128),
                         each =length(f1_list[[1]]))
                    
                  )
)

P1 <- ggplot( df1, aes(x=x, y=func))+
  geom_line(size=1.2)+
  ylab("Scale")+
  facet_wrap(.~scale, ncol=1,strip.position = "left")+
  xlab("")+
  theme_bw()+
  ggtitle("Function encoded ")+
  theme(
    
    
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    strip.placement = "outside")


df1b <- data.frame(y= do.call(c,wd1_list),
                   x= do.call(c,xwd1),
                   scale =  factor(  do.call(c, ywd1)))
df1b$color= ifelse(abs(df1b$y)>0.01, "#377eb8","#e41a1c")


P1b <- ggplot( df1b, aes(x=x, y=y, colour=color))+
  geom_point(size=2.5)+
  geom_hline(yintercept = 0)+
  facet_wrap(.~scale,scales = "free", ncol=1)+
  scale_color_manual(values= c( "#377eb8","#e41a1c") )+
  
  ggtitle("Wavelet coefficients ")+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(legend.position = "none",strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())







pi0_1 <- c( 0.99,0.99,.9, .8,.7,0.5,0.8)

df_prior1 <- data.frame(scale  = rep(c(1:log2(128)),2) ,
                        color= rep(c("#377eb8","#e41a1c"), each= log2(128)) ,
                        pi_0= c( pi0_1, 1-pi0_1),
                        sigma= rep(c(1, 0.1), each=log2(128))
)
tl <- list()
h <- 0

n_rep <- 10^4

x <- seq(-3,3, length.out=n_rep)


for ( i in 1:nrow(df_prior1)){
  
  
  tl[[ i]] <-  data.frame ( x=x,
                            y=df_prior1$pi_0[i]* dnorm(x,
                                                       sd=df_prior1$sigma[i]),
                            color=  rep(df_prior1$colo[i],n_rep),
                            scale= rep(df_prior1$scale[i],n_rep)
  )
}


res_df <- do.call(rbind,tl)

P1_prior <- ggplot( res_df, aes(x=x,y=y, colour=color))+
  geom_line(size=1.5)+
  facet_wrap(.~scale, scales = "free", ncol=1)+
  scale_color_manual(values= c( "#377eb8","#e41a1c") )+
  ylab('')+
  xlab("")+
  ggtitle("Prior on WC")+
  theme_bw()+
  theme(legend.position = "none",strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


WD_true_effect1 <- gridExtra::grid.arrange(P01,
                                           P1,P1b,P1_prior,
                                           layout_matrix = rbind(c(1,1,1),c(2,3,4)),
                                           heights=c(1,3))

 
ggsave(WD_true_effect1 , file=paste(getwd(),"/plot/wave_example.pdf",sep=""),
       width = 21,
       height = 29.7,
       units = "cm")


#### Curves------
library(susiF.alpha)
library(susieR)
set.seed(1)
X <- N3finemapping$X
 


for ( j in 1:ncol(X)){
  X[,j] <- (X[,j]-min(X[,j]))/(0.5*max(X[,j]-min(X[,j])))
}


noisy.data=list()
pos1 <- 330
pos2 <- 350#200
for ( i in 1:nrow(X))
{
  
  noise <- rnorm(length(f1  ), sd=  2)
  noisy.data [[i]] <-  X[i,pos1]*f1  +X[i,pos2]*f2   + noise
  
}

X[,pos1+1]<- X[,pos1 ]+ rnorm(nrow(X),sd=0.001)

X[,pos1-1]<- X[,pos1]+ rnorm(nrow(X),sd=0.001)
X[,pos1-2]<- X[,pos1 ]+ rnorm(nrow(X),sd=0.001)


X[,pos2+1]<- X[,pos2 ]+ rnorm(nrow(X),sd=0.0001)
X[,pos2-1]<- X[,pos2 ]+ rnorm(nrow(X),sd=0.0001)
#X[,pos2-2]<- X[,pos2 ]+ rnorm(nrow(X),sd=0.0001)

noisy.data <- do.call(rbind, noisy.data)

Y <- noisy.data
out <- susiF(Y[1:300,],X[1:300,1:500],L=2)
set.seed(1)

#### Function 1------
f1_noisy  <-  f1  + rnorm(128)
f1_noisy_list <- list()
f1_est <- out$fitted_func[[1]]
f1_est_list <- list()
wd1_list <- list()
wd1_est  <- list()
xwd1 <-list()
ywd1 <-list()
for (s in 1:log2(length( f1_noisy))){
  tt <- rep( 0, 128)
  
  
  tt_wd <- wd(tt)
  tt_wd$D[ idx_lst[[s]]]  <- wd( f1_noisy)$D[ idx_lst[[s]]]
  
  wd1_list[ idx_lst[[s]]]  <- wd( f1_noisy)$D[ idx_lst[[s]]]
  xwd1[[s]] <- (idx_lst[[s]]-mean(idx_lst[[s]]))/length(idx_lst[[s]])
  ywd1[[s]] <- rep( s,length(idx_lst[[s]]))
  f1_noisy_list[[s]] <- wr(tt_wd)
  
  
  tt <- rep( 0, 128)
  
  tt_wd$D[ idx_lst[[s]]]  <- wd( f1_est)$D[ idx_lst[[s]]]
  
  wd1_est[ idx_lst[[s]]]  <- wd( f1_est)$D[ idx_lst[[s]]]
  
  tt_wd <- wd(tt)
  tt_wd$D[ idx_lst[[s]]]  <- wd( f1_est )$D[ idx_lst[[s]]]
  f1_est_list[[s]] <- wr(tt_wd)
  
  
}


df01_noisy <- data.frame(func = c(  f1_noisy ),
                         x=  ( 1:128 )
                         
)



df_plot <- do.call( rbind, lapply(1:nrow(Y), function( i) cbind (Y[i,], 1:length(Y[i,]),rep(i, length(Y[i,]))) ))
df_plot <- data.frame(df_plot)
colnames(df_plot)<- c("y","x","ind")

set.seed(1)
noisy.data2=list()

x_1 <- sample( c(0,1,2), size=500, replace=TRUE)
for ( i in 1:500)
{
  
  noise <- rnorm(length(f1  ), sd=  0.8)
  noisy.data2 [[i]] <-  x_1[i]*f1 +   noise
  
}
noisy.data2 <- do.call(rbind, noisy.data2)

Y2 <- noisy.data2

df_plot2 <- do.call( rbind, lapply(1:nrow(Y2), function( i) cbind (Y2[i,], 1:length(Y2[i,]),rep(i, length(Y2[i,]))) ))
df_plot2 <- data.frame(df_plot2)
colnames(df_plot)<- c("y","x","ind")
p <- ggplot()

for ( i in 1: nrow (Y2)){
  dat  <- data.frame(y = Y2[i,], x = 1:length(Y2[i,]))
  p <- p + geom_line(data = dat , mapping = aes(x = x, y = y) ,alpha=.1)
}

colors <- c("SNP=1" = "palegreen3", "SNP=0" = "palegreen1", "SNP=2" = "palegreen4")

dat  <- data.frame(y =  f1, x = 1:length(Y2[i,]))
p <- p + geom_line(data = dat , mapping = aes(x = x, y = y,
                                              color = "SNP=1") ,lwd=1.5 )+
  geom_line(data = dat , mapping = aes(x = x, y = 0 *y,
                                       
                                       color = "SNP=0") ,lwd=1.5 )+
  geom_line(data = dat , mapping = aes(x = x, y = 2*y,
                                       
                                       color = "SNP=2") ,lwd=1.5 ) +
  scale_color_manual(values = colors)
P_obs2 <-p+
  theme_bw()+
  labs(color=NULL)+
  ggtitle("Observed methylation profile")+
  ylab("Methylation level")+xlab('base pair')





P01_noisy <- ggplot( df01_noisy, aes(x=x, y=func))+
  geom_line(size=2,col="green4")+
  theme_bw()+
  ylab("")+
  ggtitle("Noisy estimate")+
  theme( )






df1_noisy <- data.frame(func = c(  do.call(c,  f1_noisy_list)),
                        x= rep( 1:128,  length( f1_noisy_list) ),
                        scale =  factor(
                          rep( 1:log2(128),
                               each =length( f1_noisy_list[[1]]))
                          
                        )
)





df1_noisyb <- data.frame(y= do.call(c,wd1_list),
                         x= do.call(c,xwd1),
                         scale =  factor(  do.call(c, ywd1)))
df1_noisyb$color= ifelse(abs(df1_noisyb$y)>0.01, "#377eb8","#e41a1c")





df1_est <- data.frame(func = c(  do.call(c,  f1_est_list)),
                      x= rep( 1:128,  length( f1_est_list) ),
                      scale =  factor(
                        rep( 1:log2(128),
                             each =length( f1_est_list[[1]]))
                        
                      )
)
df1_estb <- data.frame(y= do.call(c,wd1_est),
                       x= do.call(c,xwd1),
                       scale =  factor(  do.call(c, ywd1)))
df1_estb$color= ifelse(abs(df1_estb$y)>0.01, "#377eb8","#e41a1c")






P1b_noisy <- ggplot( df1_noisyb, aes(x=x, y=y ))+
  geom_point(size=2.5,col= "#377eb8")+
  geom_hline(yintercept = 0)+
  facet_wrap(.~scale,scales = "free", ncol=1, strip.position = "left")+
  
  
  ggtitle("Noisy wavelet coefficients estimates")+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(legend.position = "none", axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.placement = "outside")






effect <- data.frame( mid= do.call(c, out$fitted_func ),
                      low= do.call(c, lapply( 1:length(out$cs), function(i)
                        out$cred_band[[i]][2,])),
                      up= do.call(c, lapply( 1:length(out$cs), function(i)
                        out$cred_band[[i]][1,])),
                      col = factor( rep( 1:length(out$cs),each=length(out$fitted_func [[1]]) )
                                    
                      ),
                      pos= rep(1:length(out$fitted_func [[1]]) , length(out$cs))
)






pi0_1est <- c( 0.96,0.99,.9, .85,.6,0.7,0.85)

df_prior1est <- data.frame(scale  = rep(c(1:log2(128)),2) ,
                           color= rep(c("#377eb8","#e41a1c"), each= log2(128)) ,
                           pi_0= c( pi0_1est, 1-pi0_1est),
                           sigma= rep(c(1, 0.1), each=log2(128))
)
tl <- list()
h <- 0

n_rep <- 10^4

x <- seq(-3,3, length.out=n_rep)


for ( i in 1:nrow(df_prior1)){
  
  
  tl[[ i]] <-  data.frame ( x=x,
                            y=df_prior1est$pi_0[i]* dnorm(x,
                                                          sd=df_prior1$sigma[i]),
                            color=  rep(df_prior1$colo[i],n_rep),
                            scale= rep(df_prior1$scale[i],n_rep)
  )
}


res_df <- do.call(rbind,tl)

P1_priorest <- ggplot( res_df, aes(x=x,y=y, colour=color))+
  geom_line(size=1.5)+
  facet_wrap( scale~., scales = "free", ncol=1, strip.position = "left")+
  scale_color_manual(values= c( "#377eb8","#e41a1c") )+
  ylab("" )+
  xlab("")+
  ggtitle("Estimated scale depend prior ")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.placement = "outside")

t_effect <-effect[which(effect$col==2),]

P1b_est <- ggplot( df1_estb, aes(x=x, y=y ))+
  geom_point(size=2.5,col= "#377eb8")+
  geom_hline(yintercept = 0)+
  facet_wrap(.~scale,scales = "free", ncol=1, strip.position = "left")+
  
  
  ggtitle("Shrunk wavelet coefficients ")+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(legend.position = "none", axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.placement = "outside")




MAP_est <- ggplot( t_effect  )+
  geom_line(aes(x = pos, y = mid,color =col), size=1.4)+
  geom_ribbon(aes(x = pos,
                  y = mid,
                  ymin = low,
                  ymax = up,
                  fill = col),
              alpha = 0.2)+
  ggtitle(" Estimated Posterior Effect") +
  xlab("base pair")+
  ylab("Effect size")+scale_color_manual(values=c("darkgreen"))+
  scale_fill_manual( values = c(   "darkgreen" ))+
  theme_bw()+
  theme(legend.position = "none" )



P_obs2
P01_noisy


P1b_noisy

P1_priorest


P1b_est


MAP_est

grid.arrange(P01_noisy, P1b_noisy)

 