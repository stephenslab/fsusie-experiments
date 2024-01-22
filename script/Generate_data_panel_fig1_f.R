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


WD_true_effect1

#ggsave(WD_true_effect1 , file="combined_plot/WD_true_effect1.pdf")













#### Function 2------
set.seed(1)
idx_lst <- gen_wavelet_indx(log2(128))


f2_list <- list()
wd2_list <- list()

xwd2 <-list()
ywd2 <-list()
for (s in 1:log2(length(f2))){
  tt <- rep( 0, 128)


  tt_wd <- wd(tt)
  tt_wd$D[ idx_lst[[s]]]  <- wd(f2)$D[ idx_lst[[s]]]
  wd2_list[ idx_lst[[s]]]  <- wd(f2)$D[ idx_lst[[s]]]
  xwd2[[s]] <- (idx_lst[[s]]-mean(idx_lst[[s]]))/length(idx_lst[[s]])
  ywd2[[s]] <- rep( s,length(idx_lst[[s]]))
  f2_list[[s]] <- wr(tt_wd)
}




df20 <-data.frame(func = c( f2),
                  x= rep( 1:128)
)
P02  <- ggplot( df20, aes(x=x, y=func))+
  geom_line(size=2, col="#6A3D9A")+
  theme_bw()+
  ylab("")+
  ggtitle("Effect to be estimated")+
  theme( )


df2 <-data.frame(func = c(do.call(c, f2_list)),
                 x= rep( 1:128,  length(f2_list) ),
                 scale =  factor(
                   rep( 1:log2(128),
                        each =length(f2_list[[1]]))
                 )
)

P2 <- ggplot( df2, aes(x=x, y=func))+
  geom_line(size=1.2)+
  ylab("Scale")+
  facet_wrap(.~scale, ncol=1,strip.position = "left")+
  xlab("")+
  theme_bw()+
  ggtitle("Function encoded")+
  theme(


    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    strip.placement = "outside")

df2b <- data.frame(y= do.call(c,wd2_list),
                   x= do.call(c,xwd2),
                   scale =  factor(  do.call(c, ywd2)))
df2b[2,1] <-0.11
df2b[1,1] <-0.11
df2b$color= ifelse(abs(df2b$y)>0.1, "#377eb8","#e41a1c")

df2b$y= ifelse(abs(df2b$y)>0.1, df2b$y,0)


P2b<-  ggplot( df2b, aes(x=x, y=y, colour=color))+
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






pi0_2 <- c( 0.99,0.85,.5, .5,.6,0.7,0.8)

df_prior2 <- data.frame(scale  = rep(c(1:log2(128)),2) ,
                        color= rep(c("#377eb8","#e41a1c"), each= log2(128)) ,
                        pi_0= c( pi0_2, 1-pi0_2),
                        sigma= rep(c(1, 0.1), each=log2(128))
)
tl <- list()
h <- 0

n_rep <- 10^4

x <- seq(-3,3, length.out=n_rep)


for ( i in 1:nrow(df_prior1)){


  tl[[ i]] <-  data.frame ( x=x,
                            y=df_prior2$pi_0[i]* dnorm(x,
                                                       sd=df_prior1$sigma[i]),
                            color=  rep(df_prior1$colo[i],n_rep),
                            scale= rep(df_prior1$scale[i],n_rep)
  )
}


res_df <- do.call(rbind,tl)

P2_prior <- ggplot( res_df, aes(x=x,y=y, colour=color))+
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



WD_true_effect2 <- gridExtra::grid.arrange(P02,
                        P2,P2b,P2_prior,
                        layout_matrix = rbind(c(1,1,1),c(2,3,4)),
                        heights=c(1,3))

WD_true_effect2

#ggsave(WD_true_effect2 , file="combined_plot/WD_true_effect2.pdf")


#### Curves------
library(susiF.alpha)
library(susieR)
set.seed(1)
X <- N3finemapping$X
for ( j in 1:ncol(X)){
  X[,j] <- (X[,j]-min(X[,j]))/(0.5*max(X[,j]-min(X[,j])))
}


noisy.data=list()
pos1 <- 100
pos2 <- 200#200
for ( i in 1:nrow(X))
{

  noise <- rnorm(length(f1  ), sd=  1)
  noisy.data [[i]] <-  X[i,pos1]*f1  +X[i,pos2]*f2   + noise

}
noisy.data <- do.call(rbind, noisy.data)

Y <- noisy.data
 dim(Y)
df_plot <- do.call( rbind, lapply(1:nrow(Y), function( i) cbind (Y[i,], 1:length(Y[i,]),rep(i, length(Y[i,]))) ))
df_plot <- data.frame(df_plot)
colnames(df_plot)<- c("y","x","ind")
p <- ggplot()

for ( i in 1: nrow (Y)){
  dat  <- data.frame(y = Y[i,], x = 1:length(Y[i,]))
  p <- p + geom_line(data = dat , mapping = aes(x = x, y = y) ,alpha=.1)
}
#dat  <- data.frame(y =  f1, x = 1:length(Y[i,]))
#p <- p + geom_line(data = dat , mapping = aes(x = x, y = y) ,lwd=1.5, col="darkgreen")
#dat  <- data.frame(y =   f2, x = 1:length(Y[i,]))
#p <- p + geom_line(data = dat , mapping = aes(x = x, y = y) ,lwd=1.5, col="#6A3D9A")

P_obs <-p+
  theme_bw()+
  ggtitle("Observed methylation profile")+
  ylab("Methylation level")+xlab('base pair')
P_obs

library(latex2exp)
#### manhattan plot ----
set.seed(12)
y <- rnorm(nrow(X),sd=4)+0.5*X[,950]+0.5*X[,920]+0.5*X[,970]
pv_list <-c()
X_tt <- X
X   [,200:400] <- X[,800:1000]
for (j in 1:ncol(X)){

  fit <-  lm(y~X[,j])
  pv_list <- c( pv_list,summary(fit)$coefficients[2,4])

}
lpv = (-log10(pv_list))

lpv[ -c(200:400) ] <- (lpv[ -c(200:400) ]/2 )^1.5


df_pv <- data.frame( lpv=lpv,
                    x= 1:length(pv_list))

df_point <- data.frame(x= c(336,360),
                       y= lpv[c(336,360)],
                       color= c("darkgreen","#6A3D9A" ),
                       label= c("Causal SNP 1","Causal SNP 2"))

pv_plot <- ggplot(df_pv)+
  geom_point(size=2,aes(y= lpv,
                       x= x,
                       color=-lpv)) +
  #geom_point(df_point, mapping=aes(x=x,y=y),colour=df_point$color ,size=4)+
  # geom_text_repel(df_point ,mapping=aes(x=x,y=y,label =factor( label)),
  #                 box.padding   = 0.35,
  #                 point.padding = 0.5 ,
  #                 nudge_x = 50,
  #                 nudge_y = -0.15,
  #                 colour=df_point$color
  #)+
  #xlim(c(0,500))+
  theme_bw()+
  theme(legend.position = "none",strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ggtitle("Marginal associations")+
  ylab(TeX("-$\\log_{10}$ p-value"))+
  xlab('SNPs')+
  xlim(c( 1,500))+
  scale_color_gradientn(colours = rainbow(3))
pv_plot
library(gridExtra)
Obs_marg <- grid.arrange(P_obs,pv_plot)
Obs_marg

#ggsave(Obs_marg , file="combined_plot/Obs_marg.pdf")


### fsusie ----
library(susiF.alpha)
library(susieR)
set.seed(1)

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
#plot_susiF(out)
pip_df <- data.frame(y =  out$pip,
                     x=1: length(out$pip),
                     color= rep("black",
                                length(out$pip)
                     )
)
pip_df$color[out$cs[[1]]]<- "darkgreen"

pip_df$color[out$cs[[2]]]<-"#6A3D9A"


x = do.call (c, lapply(1:length(out$cs), function(l)
  out$cs[[l]][which.max(out$pip[ out$cs[[l]]])]

))
y= do.call (c, lapply(1:length(out$cs), function(l)
  out$pip[out$cs[[l]][which.max(out$pip[ out$cs[[l]]])]]

))
df_label <- data.frame(x= x,
                       y=y,
                       label=c("CS1","CS2"),
                       color=factor(c("darkgreen","#6A3D9A"))
)

pip_plot <- ggplot(pip_df ,aes(x=x,y=y, colour=as.factor(color)))+
  geom_point(size=2.5)+
  scale_color_manual("Credible set",values = c( "#6A3D9A","black", "darkgreen"))+
  theme_bw()+theme(legend.position = "none",strip.background = element_blank(),
                   strip.text.x = element_blank(),
                   axis.text.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank())+
  geom_text_repel(df_label ,mapping=aes(x=x,y=y,label =factor( label)),
                   box.padding   = 0.35,
                   point.padding = 0.5 ,
                   nudge_x = 20
                   )+
  ylab("PIP")+
  xlab("SNPs")




mat_effect <- do.call(rbind, out$fitted_func )
saveRDS(mat_effect, file = "data/fig_1_data/mat_effect_orginal_space.rds")
