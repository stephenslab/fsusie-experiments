set.seed(1)

library(ggplot2)
#### geno -----


library(susieR)
data(N3finemapping)
attach(N3finemapping)

X <- N3finemapping$X

X <- X[1:20, 1:100]

library(susieR)
df <- expand.grid(x = 1:nrow(X), y =1:ncol(X))
set.seed(1)



df$z <- 0*sample(c(0,1,2), size= nrow(df), replace=TRUE)
k=0
for (  j in 1:ncol(X)){
  
  
  df$z[which(df$y==j)] <-  X[ ,j]-min(X[,j])
  if( mean( df$z[which(df$y==j)] )>0.75){
    df$z[which(df$y==j)] <- 2- df$z[which(df$y==j)]
  }
}

# default is compatible with geom_tile()
geno_plot <-  ggplot( (df), aes(x=y, y=x, fill = z)) +
  geom_raster()+
  theme_minimal()+
  scale_fill_gradient2(midpoint = 1,
                       low = "steelblue4" ,
                       mid =,"palegreen",
                       high = "#FDBF6F", space = "Lab" )+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(), legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())+
  ggtitle("Observed loci")+
  ylab("Individual")+xlab('SNP')






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


i=11
p <- ggplot()
dat  <- data.frame(y = Y[i,], x = 1:length(Y[i,]))
p <- p + geom_line(data = dat , mapping = aes(x = x, y = y) ,alpha= 1, linewidth=1.5)

#dat  <- data.frame(y =  f1, x = 1:length(Y[i,]))
#p <- p + geom_line(data = dat , mapping = aes(x = x, y = y) ,lwd=1.5, col="darkgreen")
#dat  <- data.frame(y =   f2, x = 1:length(Y[i,]))
#p <- p + geom_line(data = dat , mapping = aes(x = x, y = y) ,lwd=1.5, col="#6A3D9A")

P_obs <-p+
  theme_bw()+
  ggtitle("Observed methylation profile")+
  ylab("Methylation level")+xlab('base pair')
P_obs




f1_list <- list()
wd1_list <- list()

xwd1 <-list()
ywd1 <-list()

Y_map <- dat$y
idx_lst <-  susiF.alpha:::gen_wavelet_indx(log2(length(Y_map)))
f1 <-  dat$y
library(wavethresh)
library(ashr)
for (s in 1:log2(length(Y_map))){
  tt <- rep( 0, length(Y_map) )
  
  
  tt_wd <- wd(tt)
  tt_wd$D[ idx_lst[[s]]]  <- wd(f1)$D[ idx_lst[[s]]]
  
  wd1_list[ idx_lst[[s]]]  <- wd(f1)$D[ idx_lst[[s]]]
  xwd1[[s]] <- (idx_lst[[s]]-mean(idx_lst[[s]]))/length(idx_lst[[s]])
  ywd1[[s]] <- rep( s,length(idx_lst[[s]]))
  f1_list[[s]] <- wr(tt_wd)
  
  #  ash( wd(f1)$D[ idx_lst[[s]]])
  
  
}


df1b <- data.frame(y= do.call(c,wd1_list),
                   x= do.call(c,xwd1),
                   scale =  factor(  do.call(c, ywd1)))
df1b$color= ifelse(abs(df1b$y)>0.01, "#377eb8","#e41a1c")


P1b <- ggplot( df1b, aes(x=x, y=y, colour=color))+
  geom_point(size=2.5)+
  geom_hline(yintercept = 0)+
  facet_wrap(.~scale, ncol=1,strip.position = "left")+
  scale_color_manual(values= c( "#377eb8","#e41a1c") )+
  
  ggtitle("Observed wavelet coefficients ")+
  ylab("Scale")+
  theme_bw()+  
  theme(
    legend.position = "none",
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    strip.placement = "outside") +
  xlab("")

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
plot_susiF(out)
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






effect <- data.frame( mid= do.call(c, out$fitted_func ),
                      low= do.call(c, lapply( 1:length(out$cs), function(i)
                        out$cred_band[[i]][2,])),
                      up= do.call(c, lapply( 1:length(out$cs), function(i)
                        out$cred_band[[i]][1,])),
                      col = factor( rep( 1:length(out$cs),each=length(out$fitted_func [[1]]) )
                                    
                      ),
                      pos= rep(1:length(out$fitted_func [[1]]) , length(out$cs))
)



P_effect <- ggplot( effect  )+
  geom_line(aes(x = pos, y = mid,color = col), size=1.4)+
  geom_ribbon(aes(x = pos,
                  y = mid,
                  ymin = low,
                  ymax = up,
                  fill = col),
              alpha = 0.2)+
  ggtitle("Estimated Effect") +
  xlab("base pair")+
  ylab("Effect size")+
  geom_hline(yintercept = 0)+
  facet_wrap(.~col, ncol=1)+
  theme_bw()+
  scale_color_manual("Credible set",values = c(   "darkgreen", "#6A3D9A"))+
  
  scale_fill_manual("Credible set",values = c(   "darkgreen", "#6A3D9A"))+
  theme( strip.background = element_blank(),
         strip.text.x = element_blank(),
         legend.position="bottom")




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
 
df1b <- data.frame(y= do.call(c,wd1_list),
                   x= do.call(c,xwd1),
                   scale =  factor(  do.call(c, ywd1)))
df1b$color= ifelse(abs(df1b$y)>0.01, "#377eb8","#e41a1c")


P1 <- ggplot( df1b, aes(x=x, y=y, colour=color))+
  geom_point(size=2.5)+
  geom_hline(yintercept = 0)+
  facet_wrap(.~scale,scales = "free", ncol=1)+
  scale_color_manual(values= c( "#377eb8","#e41a1c") )+
  
  ggtitle("Wavelet coefficients ")+
  ylab("Scale")+
  theme_bw()+ 
  theme(legend.position = "none",strip.background = element_blank())+
  facet_wrap(.~scale, ncol=1,strip.position = "left")+
  xlab("")







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



 

df2 <-data.frame(func = c(do.call(c, f2_list)),
                 x= rep( 1:128,  length(f2_list) ),
                 scale =  factor(
                   rep( 1:log2(128),
                        each =length(f2_list[[1]]))
                 )
)
 
df2b <- data.frame(y= do.call(c,wd2_list),
                   x= do.call(c,xwd2),
                   scale =  factor(  do.call(c, ywd2)))
df2b[2,1] <-0.11
df2b[1,1] <-0.11
df2b$color= ifelse(abs(df2b$y)>0.1, "#377eb8","#e41a1c")

df2b$y= ifelse(abs(df2b$y)>0.1, df2b$y,0)


P2 <-  ggplot( df2b, aes(x=x, y=y, colour=color))+
  geom_point(size=2.5)+
  geom_hline(yintercept = 0)+
  facet_wrap(.~scale,scales = "free", ncol=1)+
  scale_color_manual(values= c( "#377eb8","#e41a1c") )+
  
  ggtitle("Wavelet coefficients ")+
  ylab("Scale")+
  theme_bw()+ 
  theme(legend.position = "none",strip.background = element_blank())+
  facet_wrap(.~scale, ncol=1,strip.position = "left")+
  xlab("")


P_effect
pip_plot
geno_plot
P_obs
P1b
 
ggsave(P_effect , file="data/fig_2_data/P_effect.pdf",
       width=29.7,
       height = 21,
       unit="cm")
ggsave(pip_plot  , file="data/fig_2_data/pip_plot.pdf",
       width=29.7,
       height = 21,
       unit="cm")
ggsave(geno_plot , file="data/fig_2_data/geno_plot.pdf",
       width=29.7,
       height = 21,
       unit="cm")
ggsave(P_obs , file="data/fig_2_data/observed_DNAm.pdf",
       width=29.7,
       height = 21,
       unit="cm")

ggsave(P1b , file="data/fig_2_data/raw_wc.pdf" ,
       width=29.7,
       height = 21,
       unit="cm")

mat <- do.call( rbind, out$est_pi[[1]])
row.names(mat) <- paste("scale", 0:(nrow(mat)-1))
saveRDS(mat, file = "data/fig_2_data/fitted_weight_data.rds")



