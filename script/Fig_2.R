set.seed(1)

library(ggplot2)
#### geno -----


library(susieR)
data(N3finemapping)
attach(N3finemapping)

X <- N3finemapping$X

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

#ggsave(geno_plot, file ="D:/Document/Serieux/Travail/Data_analysis_and_papers/susiF_data_simu/fig2/geno.png" )






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
#ggsave(P_obs, file ="D:/Document/Serieux/Travail/Data_analysis_and_papers/susiF_data_simu/fig2/curves.png" )



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

library(gridExtra)
Fsusie_res_plot <- grid.arrange(  P_effect ,pip_plot, ncol=1)
Fsusie_res_plot
#ggsave(Fsusie_res_plot, file ="D:/Document/Serieux/Travail/Data_analysis_and_papers/susiF_data_simu/fig2/Fsusie.png")







library(ggplot2)
library(latex2exp)
df
df= data.frame(x= 0, y=1, lab= "Functional regression model")
functional_reg_lab <- ggplot()+
  geom_text(df, mapping= aes(x=x,y=y,label = lab ))+
  geom_text(aes(x = 0, y = 0.95, label = TeX("$y= f(x)+\\epsilon$", output = "character")),
            parse = TRUE)+
  ylim(c(0,1))+



  geom_text(aes(x = 0, y = 0.5 , label = TeX("fSuSiE", output = "character")),
            parse = TRUE)+
 # geom_segment(aes(x = 0 , y = 0.9, xend = 0 , yend = 0.55),
#               arrow = arrow(length = unit(0.5, "cm")),size = 1.5)+
  theme_void()
df= data.frame(x= 0, y=1, lab= "1) Functional effects \nand credible bands")
Credban_lab <-  ggplot()+
  geom_text(df, mapping= aes(x=x,y=y,label = lab ))+

  theme_void()
df= data.frame(x= 0, y=1, lab= "2) Credible sets")

CS_lab <-  ggplot()+
  geom_text(df, mapping= aes(x=x,y=y,label = lab ))+

  theme_void()


#geno_plot
#P_obs
#P_effect
#pip_plot

#functional_reg_lab
#fSuSiE_lab_plot
#Credban_lab
#CS_lab

library(ggplot2)
library(gridExtra)
library(grid)

b = nullGrob() # per @baptiste's comment, use nullGrob() instead of rectGrob()

# grid.bezier with a few hard-coded settings
mygb = function(x,y) {
  grid.bezier(x=x, y=y, gp=gpar(fill="black"),
              arrow=arrow(type="closed", length=unit(2,"mm")))
}




png(filename="Fig2.png", units="mm", width = 400, height =200,res=300 )

grid.arrange(arrangeGrob(geno_plot,
                         P_obs, heights=c(0.5,0.5)),

             functional_reg_lab,

             arrangeGrob(Credban_lab,
                         CS_lab, heights=c(0.5,0.5)),
             arrangeGrob(P_effect,
                         pip_plot, heights=c(0.5,0.5)),
             ncol=4, widths=c(0.35, 0.2,  0.1,0.35))
vp = viewport(x = 0, y=0, width=1, height=1)
pushViewport(vp)
popViewport()

# Add top set of arrows
mygb(x=c(0.35,0.37, 0.395, 0.425), y=c(0.75,0.725,0.51,0.507))

mygb(x=c(0.35,0.37,0.395, 0.425), y=c(0.25,0.275,0.49,0.492))


mygb(x=c(0.45,0.45, 0.45, 0.45), y=(c(0.89,0.75,.6,0.51)))

mygb(x=c(0.470 ,0.505, 0.525 ,0.555), y=rev(c(0.25,0.26,0.45,0.495)))
mygb(x=c(0.470 ,0.505, 0.525 ,0.555), y=rev(c(0.75,0.74,0.55,0.505)))


# Up to "main" viewport (the "full" canvas of the main layout)
dev.off()
