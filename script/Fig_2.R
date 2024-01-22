effect <-   1.2*cos( (1:128)/128 * 3*pi )
effect[which(effect<0)]<- 0
effect[1:40]<- 0

plot( effect)
library(susiF.alpha)
library(susieR)
library(wavethresh)
set.seed(2)
data(N3finemapping)
X <- N3finemapping$X[1:100,]

set.seed(1)

obs <- list()

true_pos <- 440


for (i in 1:100 ){
  obs[[i]] <- X[i,true_pos]*effect+ rnorm(length(effect))
}

y <- do.call(rbind, obs)

table (X[,true_pos])


susie_est <- list()

for (i in which( effect>0)){
  
  
  susie_est [[i]]<- susie( X=X,y=y[,i], L=1)$sets
}

susie_est
#88
#91
#94
#100
library(ggplot2)
#88
df_effect <-data.frame(x=1:128,
                       y=effect)
P00 <- ggplot( df_effect, aes(x=x, y=y))+
  geom_point()+
  
  ylim(c(-0.1,1.3))+
  theme_classic()+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1.2),axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank())
P01 <- P00+geom_point(x=88,y=effect[88], col="red", shape=21,size=3)+
  
  geom_hline(yintercept = 0)+
  xlab("CpG")+
  ylab("Effect ") 

tt <- susie( X=X,y=y[,88], L=1)

df_pip1 <-data.frame(x=1:length(tt$pip),
                     y=tt$pip)
df2   <-data.frame(x= tt$sets$cs$L1,
                   y=tt$pip[tt$sets$cs$L1])
df3   <-data.frame(x= true_pos,y=tt$pip[true_pos])

P11 <- ggplot( )+
  geom_point(df_pip1, mapping=aes(x=x, y=y))+
  xlab("SNP index")+
  ylab("PIP")+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_point(df2 ,  mapping=aes(x=x,y=y), col="lightblue3",size=3)+
  geom_point(df3 ,  mapping=aes(x=x,y=y), col="red", shape=21,size=3)

P11




P02 <- P00+geom_point(x=91,y=effect[91], col="red", shape=21,size=3)+
  ylab("")+
  xlab("CpG")+
  ylab("Effect ")+
  geom_hline(yintercept = 0)+
  theme( axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank())
tt <- susie( X=X,y=y[,91], L=1)

df_pip2 <-data.frame(x=1:length(tt$pip),
                     y=tt$pip)
df2   <-data.frame(x= tt$sets$cs$L1,
                   y=tt$pip[tt$sets$cs$L1])
df3   <-data.frame(x= true_pos,y=tt$pip[true_pos])
P21 <- ggplot( )+
  geom_point(df_pip2, mapping=aes(x=x, y=y))+
  xlab("SNP index")+
  ylab("PIP")+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  
  geom_point(df2 ,  mapping=aes(x=x,y=y), col="lightblue3",size=3)+
  geom_point(df3 ,  mapping=aes(x=x,y=y), col="red", shape=21,size=3)


P03 <- P00+geom_point(x=94,y=effect[94], col="red", shape=21,size=3)+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  
  geom_hline(yintercept = 0)+
  xlab("CpG")+
  ylab("Effect ") 

tt <- susie( X=X,y=y[,94], L=1)

df_pip3 <-data.frame(x=1:length(tt$pip),
                     y=tt$pip)
df2   <-data.frame(x= tt$sets$cs$L1,
                   y=tt$pip[tt$sets$cs$L1])
df3   <-data.frame(x= true_pos,y=tt$pip[true_pos])
P31 <-ggplot( )+
  geom_point(df_pip3, mapping=aes(x=x, y=y))+
  xlab("SNP index")+
  ylab("PIP")+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  
  geom_point(df2 ,  mapping=aes(x=x,y=y), col="lightblue3",size=3)+
  geom_point(df3 ,  mapping=aes(x=x,y=y), col="red", shape=21,size=3)

P04 <- P00+geom_point(x=100,y=effect[100], col="red", shape=21,size=3)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  
  geom_hline(yintercept = 0)+
  xlab("CpG")+
  ylab("Effect ") 

tt <- susie( X=X,y=y[,100], L=1)

df_pip4 <-data.frame(x=1:length(tt$pip),
                     y=tt$pip)
df2   <-data.frame(x= tt$sets$cs$L1,
                   y=tt$pip[tt$sets$cs$L1])
df3   <-data.frame(x= true_pos,y=tt$pip[true_pos])
P41 <-  ggplot( )+ theme_classic()+
  geom_point(df_pip4, mapping=aes(x=x, y=y))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  xlab("SNP index")+
  ylab("PIP")+
  
  
  geom_point(df2 ,  mapping=aes(x=x,y=y), col="lightblue3",size=3)+
  geom_point(df3 ,  mapping=aes(x=x,y=y), col="red", shape=21,size=3)


tt <- susiF(y, X, L=1)
tt$cs
#ten column where removed so we need to remake the plot for the explanation (other t)
tt$pip <- rep( 0,length(tt$pip))
tt$pip[true_pos] <-1

df_pip5 <-data.frame(x=1:length(tt$pip),
                     y=tt$pip)
df_est_f<- data.frame(y=c( tt$fitted_func[[1]],
                           tt$cred_band[[1]][1,],
                           tt$cred_band[[1]][2,]
),
x= rep( 1:128, 3),
type=factor(rep(1:3, each=128)))



df_est_f<- data.frame(y=c( tt$fitted_func[[1]],
                           tt$cred_band[[1]][1,],
                           tt$cred_band[[1]][2,]
),
x= rep( 1:128, 3),
type=factor(rep(1:3, each=128)))


ggplot( )+
  geom_point(df_effect,  mapping=aes(x=x, y=y))+
  geom_line( df_est_f[which(df_est_f$type==1),],
             mapping=aes(x=x, y=y,linetype="longdash"), 
             col="lightblue3",
             size=1.3)+
  geom_line( df_est_f[which(df_est_f$type==2),],
             mapping=aes(x=x, y=y,linetype="solid"),
             col="lightblue3",
             size=1.3)+
  geom_line( df_est_f[which(df_est_f$type==3),],
             mapping=aes(x=x, y=y,linetype="solid"),
             col="lightblue3",
             size=1.3)+
  
  xlab("")+
  ylab("")+
  theme_classic()

P05 <-  ggplot( )+
  geom_point(df_effect,  mapping=aes(x=x, y=y))+
  geom_line( df_est_f[which(df_est_f$type==1),],
             mapping=aes(x=x, y=y,linetype="longdash"), 
             col="lightblue3",
             size=1.3)+
  geom_line( df_est_f[which(df_est_f$type==2),],
             mapping=aes(x=x, y=y,linetype="solid"),
             col="lightblue3",
             size=1.3)+
  geom_line( df_est_f[which(df_est_f$type==3),],
             mapping=aes(x=x, y=y,linetype="solid"),
             col="lightblue3",
             size=1.3)+
  
  xlab("CpG")+
  ylab("Effect ")+
  theme_classic()+
  geom_hline(yintercept = 0)+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
P51 <- ggplot( df_pip5, aes(x=x, y=y))+
  geom_point()+
  xlab("SNP index")+
  ylab("PIP")+
  theme_classic()+
  geom_point(x= true_pos,y=tt$pip[true_pos], col="lightblue3",size=3)+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1.2),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank())+
  
  geom_point(x= true_pos,y=tt$pip[true_pos], col="red", shape=21,size=3)

library(cowplot)
library(gridExtra)
grid.arrange(#P01, P11,
             P02,P21,
              P03, P31,
              #P04,P41,
             P05,P51, ncol=2)
