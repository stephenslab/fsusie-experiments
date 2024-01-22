effect <-   0.2*cos( (1:128)/128 * 3*pi )
effect[which(effect<0)]<- 0
effect[1:40]<- 0

plot( effect)
library(susiF.alpha)
library(susieR)
library(wavethresh)
set.seed(2)
data(N3finemapping)
X <- N3finemapping$X[1:60,]

set.seed(1)

obs <- list()

true_pos <- 350

table (X[,true_pos])
X[,true_pos] <- (X[,true_pos] - min(X[,true_pos] ))
table (X[,true_pos])
for (i in 1:60 ){
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
idx0 <- idx[order(idx)]
idx  <-idx0 
plot(pos[idx],effect[50:120][idx])
 
susie_est <- list()
h=1
pips <- list()
y <- do.call(rbind, obs)
y <- y[, 50:120]
for (i in idx){
  
  
  susie_est [[h]]<- susie( X=X,y=y[,i], L=1)$sets
  pips [[h]]    <-susie( X=X,y=y[,i], L=1)$pip
  h=h+1
}

susie_est

colnames(X) <- 1:ncol(X)

susif_res <- susiF(y[,idx0] , X, L=1 ,pos=pos[idx0] )
susif_res$cs
plot(pos[idx0],effect[50:120][idx0])

lines(x=susif_res$outing_grid , y=susif_res$fitted_func[[1]])

lines(x=susif_res$outing_grid , y=susif_res$cred_band[[1]][1,])
lines(x=susif_res$outing_grid , y=susif_res$cred_band[[1]][2,])

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

 
SNP_idx <-  unique (do.call ( c, lapply( 1:length(susie_est), function (i){
  susie_est[[i]]$cs$L1
}
                                                        
  )))

SNP_idx  <- SNP_idx[order(SNP_idx )]

image( abs(cor(X[,SNP_idx])))
 


library(ggplot2)
ggplot(df, aes(x=x, y=y, col=as.factor(Genotype)))+
  geom_point(alpha=.2)+
  geom_smooth(se = FALSE)



df1 <- data.frame ( y=effect[50:120][idx0]+0.25, x= .5+pos[idx0])

df2 <- data.frame ( y=2*effect[50:120][idx0]+0.25, x= .5+pos[idx0])


df3 <- data.frame ( y=0*effect[50:120][idx0]+0.25, x= .5+pos[idx0])

pos
which ( effect[50:120]>0)
 

is.detected <- !do.call(c,lapply(1:length(susie_est), function (i){
  is.null(susie_est[[i]]$cs)
}))

df_label <- data.frame(lab= letters[1:15],
                       x= unique(df$x)[which(idx >16 & idx <57)], 
                       y=rep(-0.051, 15),
                       
                       col= ifelse(is.detected[which(idx >16 & idx <57)] ==TRUE, "lightblue4",
                                   "orange4"))


 
pos_SNP_plot <-seq(0,max ( unique(df$x)),
                   length.out=length(SNP_idx))

df_SNP_pos <- data.frame( x= pos_SNP_plot,
                          lab = factor(1:length(SNP_idx)),
                          y=rep(-0.12, length(SNP_idx))
                         )

idx_causal_SNP <- which(SNP_idx==350)


library(ggplot2)
library(ggrepel)
P0 <- ggplot()+
  geom_point(data=df, aes(x=x, y=y, col=as.factor(Genotype)),alpha=.4)+
  geom_line (data=df1 , aes ( x=x, y=y ),col="green2",size=2)+
  geom_line (data=df2 , aes ( x=x, y=y ),col="slateblue1",size=2)+
  geom_line (data=df3 , aes ( x=x, y=y ),col="tomato",size=2)+
  geom_label (data=df_label, aes (x=x,y=y, label=lab,fill=col))+
  geom_point  (data=df_SNP_pos, aes (x=x,y=y ), size=3)+
   geom_point(  aes( x=df_SNP_pos$x[idx_causal_SNP] ,
                y=df_SNP_pos$y[idx_causal_SNP] ) , 
                col="red", size=3) +
  xlab("")+ylab("")+
  scale_colour_manual(values= c( "1"= "green2",
                                 "2"="slateblue1",
                                 "0"="tomato"))+
  theme_classic()+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
P0 

tt <- susie_est [[9]] 
tt0 <- pips[[9]][SNP_idx ]



df_pip1 <-data.frame(x=pos_SNP_plot,
                     y=tt0)
df2   <-data.frame(x= pos_SNP_plot [ which ( SNP_idx %in% tt$cs$L1)] ,
                   y=tt0[which ( SNP_idx %in% tt$cs$L1)])
df3   <-data.frame(x= pos_SNP_plot[idx_causal_SNP],
                   y=tt0[idx_causal_SNP])

P11 <- ggplot( )+
  geom_point(df_pip1, mapping=aes(x=x, y=y))+
  xlab("SNP index")+
  ylab("PIP")+
  ggtitle("Fine mapping on a")+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_point(df2 ,  mapping=aes(x=x,y=y), col="lightblue3",size=3)+
  geom_point(df3 ,  mapping=aes(x=x,y=y), col="red", shape=21,size=3)

P11




 


tt <- susie_est [[11]] 
tt0 <- pips[[11]][SNP_idx ]



df_pip1 <-data.frame(x=pos_SNP_plot,
                     y=tt0)
df2   <-data.frame(x= pos_SNP_plot [ which ( SNP_idx %in% tt$cs$L1)] ,
                   y=tt0[which ( SNP_idx %in% tt$cs$L1)])
df3   <-data.frame(x= pos_SNP_plot[idx_causal_SNP],
                   y=tt0[idx_causal_SNP])

P12 <- ggplot( )+
  geom_point(df_pip1, mapping=aes(x=x, y=y))+
  xlab("SNP index")+
  
  ggtitle("Fine mapping on c")+
  ylab("PIP")+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_point(df2 ,  mapping=aes(x=x,y=y), col="lightblue3",size=3)+
  geom_point(df3 ,  mapping=aes(x=x,y=y), col="red", shape=21,size=3)
 
P12







tt <- susie_est [[13]] 
tt0 <- pips[[13]][SNP_idx ]

df_pip1 <-data.frame(x=pos_SNP_plot,
                     y=tt0)
df2   <-data.frame(x= pos_SNP_plot [ which ( SNP_idx %in% tt$cs$L1)] ,
                   y=tt0[which ( SNP_idx %in% tt$cs$L1)])
df3   <-data.frame(x= pos_SNP_plot[idx_causal_SNP],
                   y=tt0[idx_causal_SNP])

P13 <- ggplot( )+
  geom_point(df_pip1, mapping=aes(x=x, y=y))+
  xlab("SNP index")+
  ylab("PIP")+
  
  ggtitle("Fine mapping on d")+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_point(df2 ,  mapping=aes(x=x,y=y), col="lightblue3",size=3)+
  geom_point(df3 ,  mapping=aes(x=x,y=y), col="red", shape=21,size=3)

P13

diff_bot <- susif_res$cred_band[[1]][2,]- susif_res$fitted_func[[1]]
df_est_f<- data.frame(y=c(   (susif_res$fitted_func[[1]]),
                             (susif_res$cred_band[[1]][1,]-diff_bot),
                             (susif_res$cred_band[[1]][2,]+diff_bot)
                                                    ),
                      x= rep( susif_res$outing_grid, 3),
                      type=factor(rep(1:3, each=32)))
 

df_effect <- data.frame ( y=  effect[50:120][idx0], 
                    x=  pos[idx0])
 

P21 <-ggplot( )+
  geom_line(df_effect,  mapping=aes(x=x, y=y))+
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
  ggtitle("FSuSiE effect estimate")+
  xlab("")+
  ylab("")+
  theme_classic()+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
P21 





tt <- susif_res  
tt0 <- rep(0,length( SNP_idx  ))
# susif reindex the ouptut so given that 42 columns are excluded a bit annoying to reindex
tt0[9]<-1
#anyhow check in susif_res$cs you will see that full weight in the causal part
df_pip1 <-data.frame(x=pos_SNP_plot,
                     y=tt0)
df2   <-data.frame(x= pos_SNP_plot [9] ,
                   y=tt0[9])
df3   <-data.frame(x= pos_SNP_plot[idx_causal_SNP],
                   y=tt0[idx_causal_SNP])
P22 <- ggplot( )+
  geom_point(df_pip1, mapping=aes(x=x, y=y))+
  xlab("SNP index")+
  ylab("PIP")+
  
  ggtitle("Fine mapping with FSuSiE")+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_point(df2 ,  mapping=aes(x=x,y=y), col="lightblue3",size=3)+
  geom_point(df3 ,  mapping=aes(x=x,y=y), col="red", shape=21,size=3)

P22

library(gridExtra)


grid_plot <- ggdraw()+
  
  draw_plot(P0,  x = 0.01 , y = .6, width = .99, height = .39)+
  draw_plot(P11,  x = 0.01 , y = .4, width = .49, height = .2 )+
  draw_plot(P12,  x = 0.01 , y = .2, width = .49, height = .2 )+
  draw_plot(P13,  x = 0.01 , y = .0, width = .49, height = .2 )+
  
  draw_plot(P21,  x = 0.5 , y = .3, width = .49, height = .3)+
  
  draw_plot(P22,  x = 0.5 , y = .0, width = .49, height = .3)
 
grid_plot


ggsave(grid_plot , file="plot/testfig2.pdf",
       width = 29.7,
       height = 21,
       units = "cm"
)
