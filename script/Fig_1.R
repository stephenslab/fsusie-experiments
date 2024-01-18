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

true_pos <- 700
for (i in 1:100 ){
  obs[[i]] <- X[i,true_pos]*effect+ rnorm(length(effect))
}

y <- do.call(rbind, obs)




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
df3   <-data.frame(x= 700,y=tt$pip[700])

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
df3   <-data.frame(x= 700,y=tt$pip[700])
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
df3   <-data.frame(x= 700,y=tt$pip[700])
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
df3   <-data.frame(x= 700,y=tt$pip[700])
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
tt$pip[700] <-1

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
  geom_point(x= 700,y=tt$pip[700], col="lightblue3",size=3)+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1.2),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank())+
  
  geom_point(x= 700,y=tt$pip[700], col="red", shape=21,size=3)

library(cowplot)
library(gridExtra)
grid.arrange(P01, P11,
             P02,P21,
             P03, P31,
             P04,P41,
             P05,P51, ncol=2)
#### methyaoltion plot ------

path <-  getwd()


m_1182_resid_plot_box <- read.delim("D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_1_data/tad1182_raw.beta.tsv")
logit <- function( p)
{
  
  out <- log( p/(1-p))
  
  return(out)
}

logit(0.2)

library(dplyr)
library(ggplot2)
df <- rbind(m_1182_resid_plot_box%>%
              mutate(type = "mQTL"))%>%
  mutate( `rs35060601` = ifelse(X ==2, "GG", "AA" ), peak_pos = as.numeric(peak))

df$rs35060601[which(df$X==1)]<- "AG"

df$dummy_pos <-  rep(NA, nrow(df))

tt <- df$peak_pos[order(df$peak_pos)]
for ( k in 1:length(unique(tt))){
  df$dummy_pos[which(df$peak_pos==unique(tt)[k] )] <- k
  
  
}

df$genotype <- factor(df$rs35060601)
df$dummy_pos2 <-  df$dummy_pos

df$dummy_pos <-  factor(df$dummy_pos )
#df <-df[ -which( df$dummy_pos2==21),]


est <-  c()
tt <- df$peak_pos[order(df$peak_pos)]
pv <- c()
for ( k in 1:length(unique(tt))){
  fit <-  lm( df$count[which(df$peak_pos==unique(tt)[k] )]  ~  as.factor(df$X[which(df$peak_pos==unique(tt)[k] )]  ))
  fit2 <-  lm( df$count[which(df$peak_pos==unique(tt)[k] )]  ~  1 )
  
  ano <-  anova(fit,fit2)
  fit <-  lm( df$count[which(df$peak_pos==unique(tt)[k] )]  ~  as.factor(df$X[which(df$peak_pos==unique(tt)[k] )]  ))
  est <-  c(est,summary( fit)$coefficients[2,1]/summary( fit)$coefficients[2,2])
  #pv <- c(pv,summary( fit)$coefficients[2,4])
  pv <- c(pv, ano$`Pr(>F)`[2])
}


t_df <- data.frame(pv=pv, est=est,x=df$dummy_pos )

P2 <- ggplot( t_df, aes(y=est,x=x))+
  geom_line()+
  geom_point() +theme_classic( )+
  geom_hline(yintercept = 0)+
  ylab("Z-score")+
  xlab("")+theme(panel.border = element_rect(colour = "black", fill=NA, size=1.2),
                 axis.text.x=element_blank())
P2

P3 <- ggplot( t_df, aes(y=-log10(pv),x=x))+
  geom_line()+
  geom_point() +theme_classic( )+
  geom_hline(yintercept = 0)+
  ylab("-log10(pv)")+
  xlab("")+theme(panel.border = element_rect(colour = "black", fill=NA, size=1.2),
                 axis.text.x=element_blank())
P3



P1_1 <-ggplot(df[which( df$dummy_pos2> 35 & df$dummy_pos2<50),] , aes(x =  as.factor(dummy_pos2),  y = count   ))+
  geom_boxplot(aes(fill=genotype),position=position_dodge(1 ),
               #size = 0.5 ,
               alpha = 0.5,coef = 6)+
  
  
  stat_summary(fun=median,
               geom="line",
               aes(group = X,color = genotype,x  = as.factor(dummy_pos2)),
               size = 1.1,alpha=0.5)+
  theme_classic( )+
  ylab("Methylation level ")+
  xlab("CpG")+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1.2),
                     legend.position="none",
                     axis.text.x =   element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.text.y =   element_blank(),
                     axis.ticks.y=element_blank())

P1_1




P1_2 <-ggplot(df[which( df$dummy_pos2>67 & df$dummy_pos2<78),] , aes(x =  as.factor(dummy_pos2),  y = count   ))+
  geom_boxplot(aes(fill=genotype),position=position_dodge(1 ),
               #size = 0.5 ,
               alpha = 0.5,coef = 6)+
  
  ylim(c(0,.17))+
  stat_summary(fun=median,
               geom="line",
               aes(group = X,color = genotype,x = as.factor(dummy_pos2)),
               size = 1.1,alpha=0.5)+
  theme_classic( )+
  ylab("Methylation level ")+
  theme_classic()+
  xlab("CpG")+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1.2),
                     legend.position="none",
                     axis.text.x = element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.text.y =   element_blank(),
                     axis.ticks.y=element_blank())


P1_2


#sP1_1

df$alpha= rep(0.2,nrow(df))

df$alpha[c(38:47, 141:146)]<- 1
#df$rs35060601 <- df$genotype


P1_prime <- ggplot(df[which(  df$dummy_pos2>22 &  df$dummy_pos2<85),] , aes(x =  as.factor(dummy_pos2),  y = count   ))+
  stat_summary(fun=median,
               geom="line",
               aes(group = X,color = genotype,x = dummy_pos),
               size = 0.8)+
  theme_classic( )+
  theme( 
    legend.position = "bottom")

P1 <-ggplot(df[which(  df$dummy_pos2>22 &  df$dummy_pos2<85),] , aes(x =  as.factor(dummy_pos2),  y = count   ))+
  # geom_boxplot(aes(fill=genotype), alpha=0.5,position=position_dodge(1 ),
  
  #              coef = 6)+
  
  
  stat_summary(fun=median,
               geom="line",
               aes(group = X,color = genotype,x = dummy_pos),
               size = 0.8)+
  
  theme_classic( )+
  
  ylab("Methylation level")+
  xlab("CpG")+
  theme_classic()+
  # geom_segment(aes(x =17, y = 0.25, xend = 25, yend = 0.25),size=1.5)+
  # geom_text(x=21, y=0.33, label="A" )+
  # geom_segment(aes(x =  46, y = 0.25, xend =55, yend = 0.25),size=1.5)+
  # geom_text(x=50.5, y=0.33, label="B" )+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1.2),
         legend.position = "none",
         axis.text.x =   element_blank(),
         axis.ticks.x=element_blank(),
         axis.text.y =   element_blank(),
         axis.ticks.y=element_blank()
  )

legend <- cowplot::get_legend(P1_prime)

library(ggridges)
set.seed(2)

pv [69:77]<-  pv [69:77]/10
pv1 <- -log10(pv)[23:84]+2

pv1[which.max(pv1)]<-6.2
pv1[22]<-7.8
pv1[21]<-6.8
pv1[20]<-5.8
pv2 <- pv1-  runif( length(pv1), max= min(pv1)/2)
tt0 <- pv2[23]
tt1 <- pv2[22]
pv2[23] <-tt1
pv2[22]<-tt0

pv3 <- pv2- runif( length(pv1), max= max(pv2)/5)
pv4 <- pv3+ runif( length(pv1), max= max(pv3)/5)

pv5 <- pv4- runif( length(pv1), max= max(pv4)/5)

pv6 <- pv5+ runif( length(pv1), max= max(pv5)/5)

pv7 <- pv6- runif( length(pv1), max= max(pv6)/5)

pv8 <- pv7+ runif( length(pv1), max= max(pv7)/5)

pv9 <- pv8- runif( length(pv1), max= max(pv8)/5)

pv10 <- pv9+ runif( length(pv1), max= max(pv9)/5)


tt <- data.frame(pv =(c(pv1-min(pv1) ,
                        pv2-min(pv2) ,
                        pv3-min(pv3) ,
                        pv4-min(pv4) ,
                        pv5-min(pv5) ,
                        pv6-min(pv6) ,
                        pv7-min(pv7) ,
                        pv8-min(pv8) ,
                        pv9-min(pv9) ,
                        pv10-min(pv10) ) 
                      
                      
                      
),
x= rep( c(1:length( 23:84 )),10),
#  type= factor(rep(c("rs35060601", "rs1595553879", "rs1963838328"),each=length( 23:84 ))))
type= factor(rep(c("SNP1", 
                   "SNP2",
                   "SNP3",
                   "SNP4",
                   "SNP5",
                   "SNP6",
                   "SNP7",
                   "SNP8",
                   "SNP9",
                   "SNP10" ),
                 each=length( 23:84 ))))

mh_plot <- ggplot( tt, aes (x=as.factor(x),
                            y=type ,
                            height = pv,
                            group=type)
)+
  geom_ridgeline(fill="lightblue", alpha=0.5,scale=2)+
  theme_classic( )+
  ylab("-log10 pvalue")+
  xlab("")+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.2),
        axis.text.x =   element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y =   element_blank(),
        axis.ticks.y=element_blank()
        # axis.text.x = element_blank(),
        #axis.text.y = element_text(angle =45, vjust = 0, hjust= 0)
  )#+scale_y_discrete(guide = guide_axis(check.overlap = TRUE))
mh_plot

#ggdraw() +draw_plot(mh_plot, x = 0, y = 0.8, width = 1, height = 0.2)+
#  draw_plot(P1, x = 0, y = 0.4, width = 1, height = 0.4) +
#  draw_plot(P1_1, x = 0, y = 0, width = .45, height = .4) +
#  draw_plot(P1_2, x = 0.5, y = 0, width = .45, height = .4) +
#  draw_plot_label(label = c(" ", "A", "B"), size = 15,
#                  x = c(0, 0 , 0.5), y = c(1, 0.4, 0.4))






P1_all <-  
  
  
  cowplot::ggdraw() +
  coord_equal(xlim = c(0, 20), ylim = c(0, 20), expand = FALSE) +
  annotation_custom(ggplotGrob(P1), xmin = 0, xmax = 20, ymin = 10, ymax = 20) +
  
  annotation_custom(ggplotGrob(P1_1), xmin = 0, xmax = 10, ymin = 0, ymax = 10)+
  geom_segment(aes(x =  5.5, xend = 0.9775   , y = 15    , yend =9.7 ), color = "grey64", size = .666) +
  geom_segment(aes(x =  8.5, xend = 9.6775, y = 15      , yend =9.7 ), color = "grey64", size = .666)+
  geom_segment(aes(x =  8.5, xend = 5.5 , y = 15      , yend =15  ), color = "grey64", size = .666)+
  geom_segment(aes(x =  8.5, xend = 5.5 , y = 18      , yend =18  ), color = "grey64", size = .666)+
  geom_segment(aes(x =  8.5, xend = 8.5 , y = 18.03  , yend =15  ), color = "grey64", size = .666)+
  geom_segment(aes(x =  5.5, xend = 5.5 , y = 18.03  , yend =15  ), color = "grey64", size = .666)+
  
  
  
  annotation_custom(ggplotGrob(P1_2), xmin = 10, xmax = 20, ymin = 0, ymax = 10)+
  geom_segment(aes(x =  14.65  , xend = 11  , y =11.2, yend =9.7 ), color = "grey64", size = .666)+
  geom_segment(aes(x =  17.5   , xend = 19.725 , y = 11.2 , yend =9.7 ), color = "grey64", size = .666)+ 
  geom_segment(aes(x =  14.65  , xend =17.5   , y = 13   , yend =13  ), color = "grey64", size = .666)+
  geom_segment(aes(x =  17.5   , xend =17.5   , y = 11.18 , yend =13.03  ), color = "grey64", size = .666)+
  geom_segment(aes(x =  14.65  , xend =14.65  ,  y = 11.18 , yend =13.03  ), color = "grey64", size = .666)+
  geom_segment(aes(x =  14.65  , xend =17.5  , y = 11.2 , yend =11.2  ), color = "grey64", size = .666)      



Fig1 <- ggdraw() +
  draw_plot(mh_plot, x = 0, y = 0.8, width = 0.5, height = 0.2)+
  draw_plot(P1_all, x = 0, y = 0.05, width = 0.5, height = 0.75) +
  # draw_plot(P1_1, x = 0, y = 0.05, width = .25, height = .4) +
  #draw_plot(P1_2, x = 0.25, y = 0.05, width = .25, height = .4) +
  draw_plot(legend, x = 0.0, y = 0.0 , width = .5, height = .05) +
  
  draw_plot(P01, x = 0.5, y = 0.8, width = 0.25, height = 0.2)+
  
  draw_plot(P11, x = 0.75, y = 0.8, width = 0.25, height = 0.2)+
  draw_plot(P02, x = 0.5,  y = 0.6, width = 0.25, height = 0.2)+
  draw_plot(P21, x = 0.75, y = 0.6, width = 0.25, height = 0.2)+
  draw_plot(P03, x = 0.5,  y = 0.4, width = 0.25, height = 0.2)+
  draw_plot(P31, x = 0.75, y = 0.4, width = 0.25, height = 0.2)+
  draw_plot(P04, x = 0.5,  y = 0.2, width = 0.25, height = 0.2)+
  draw_plot(P41, x = 0.75, y = 0.2, width = 0.25, height = 0.2)+
  draw_plot(P05, x = 0.5,  y = 0.0, width = 0.25, height = 0.2)+
  draw_plot(P51, x = 0.75, y = 0.0, width = 0.25, height = 0.2)+
  
  draw_plot_label(label = c("A ", "B", "C",  "D" , "E",   "F",  "G"), size = 12,
                  x = c(0,     0 , 0.5 , 0.5 ,  0.5 ,  0.5 , 0.5 ), 
                  y = c(1,  0.775,  1,    0.8,   0.6,   0.4,  0.2))

library(gridExtra)
library(grid)
Fig1
library(ggpubr)
y = c(0.0081  , 0.5, 0.99  )
x = c(0.5, 0.5, 0.5   )
id = c(1, 1    ,1  )


grid.polygon(x , y , id=id)


ggsave(Fig1 , file="plot/Fig1.pdf",
       width = 29.7,
       height = 21,
       units = "cm"
)
