set.seed(2)
library(reshape2)
L <- sample(1:5, size =1)#actual number of effect
lf <-  list()
R=2^7-2
for(l in 1:L){
  
  block_cov_generate = function(corr=0.9, v=1){
    R = dim(corr)[1]
    cov = diag(sqrt(v),R) %*% corr %*% diag(sqrt(v),R)
    return(cov)
  }
  
  decay_corr_generate = function(R =2^7-2,# # sites
                                 num_block = 3,# #islands within the region
                                 off_diag_corr = 0.99,
                                 decay =0.995#
  ){
    corr_b = matrix(0,nrow=R,ncol=R)
    region_mat = matrix(off_diag_corr,R/num_block,R/num_block)
    for ( i in 1: nrow(region_mat)){
      for ( j in  1:ncol(region_mat)){
        tt <- region_mat[i,j]
        region_mat[i,j] <- tt*decay ^( abs(i-j))
        region_mat[j,i] <- tt*decay ^( abs(i-j))
        
      }
    }
    for(i in 1:num_block){
      corr_b[seq(((i-1)*(R/num_block)+1),i*(R/num_block)), seq(((i-1)*(R/num_block)+1),i*(R/num_block))] = region_mat
    }
    
    diag(corr_b) = 1
    return(corr_b)
  }
  corr <-    decay_corr_generate()
  cov = block_cov_generate(corr =   corr, v = 0.1)
  
  
  lf[[l]] <- dae::rmvnorm(mean = rep(0,R),V = cov) #functional effect for effect l
}

plot( lf[[1]])


df_plot_decay  <- data.frame(y= lf[[1]], x= 1:length(lf[[1]])) 

Pf_decay  <-  ggplot( df_plot_decay, aes( x=x,y=y))+
  geom_line()+theme_void() +  # Minimal theme
  theme(   # Rotate x axis texts
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank() ,
    legend.position = "none")
Pf_decay
L <- sample(1:5, size =1)#actual number of effect
lf <-  list()
R=(2^7) -2
for(l in 1:L){
  
  block_cov_generate = function(corr=0.9, v=1){
    R = dim(corr)[1]
    cov = diag(sqrt(v),R) %*% corr %*% diag(sqrt(v),R)
    return(cov)
  }
  
  decay_corr_generate = function(R =(2^7) -2,# # sites
                                 num_block = 3,# #islands within the region
                                 off_diag_corr = 0.99 ,
                                 decay =0.95 #
  ){
    corr_b = matrix(0,nrow=R,ncol=R)
    region_mat = matrix(off_diag_corr,R/num_block,R/num_block)
    for ( i in 1: nrow(region_mat)){
      for ( j in  1:ncol(region_mat)){
        tt <- region_mat[i,j]
        region_mat[i,j] <- tt*decay ^( abs(i-j))
        region_mat[j,i] <- tt*decay ^( abs(i-j))
        
      }
    }
    for(i in 1:num_block){
      corr_b[seq(((i-1)*(R/num_block)+1),i*(R/num_block)), seq(((i-1)*(R/num_block)+1),i*(R/num_block))] = region_mat
    }
    
    diag(corr_b) = 1
    return(corr_b)
  }
  corr <-    decay_corr_generate()
  cov = block_cov_generate(corr =   corr, v = 0.1)
  
  
  lf[[l]] <- dae::rmvnorm(mean = rep(0,R),V = cov) #functional effect for effect l
}


 image (cov)
 cor_matrix <- sqrt(corr)
 melted_cor_matrix <- melt(cor_matrix)
 
 P_decay <-  ggplot(data = melted_cor_matrix, aes(Var1, Var2, fill = value)) +
   geom_tile(color = "white") +  # Use tiles, with white borders
   scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                        midpoint = 0, limit = c(-1,1), space = "Lab", 
                        name="Pearson\nCorrelation") +
   theme_void() +  # Minimal theme
 
   # Minimal theme
  theme(   # Rotate x axis texts
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank() ,
    legend.position = "none") +
  coord_fixed()  # Ensure the tiles are square
 P_decay
 
 
 set.seed(4) 
 
 L <- sample(1:5, size =1)#actual number of effect
 lf <-  list()
 R=(2^7) -2
 for(l in 1:L){
   
   block_cov_generate = function(corr=0.9, v=1){
     R = dim(corr)[1]
     cov = diag(sqrt(v),R) %*% corr %*% diag(sqrt(v),R)
     return(cov)
   }
   
   decay_corr_generate = function(R =(2^7) -2,# # sites
                                  num_block = 3,# #islands within the region
                                  off_diag_corr =1,
                                  decay =1#
   ){
     corr_b = matrix(0,nrow=R,ncol=R)
     region_mat = matrix(off_diag_corr,R/num_block,R/num_block)
     for ( i in 1: nrow(region_mat)){
       for ( j in  1:ncol(region_mat)){
         tt <- region_mat[i,j]
         region_mat[i,j] <- tt*decay ^( abs(i-j))
         region_mat[j,i] <- tt*decay ^( abs(i-j))
         
       }
     }
     for(i in 1:num_block){
       corr_b[seq(((i-1)*(R/num_block)+1),i*(R/num_block)), seq(((i-1)*(R/num_block)+1),i*(R/num_block))] = region_mat
     }
     
     diag(corr_b) = 1
     return(corr_b)
   }
   corr <-    decay_corr_generate()
   cov = block_cov_generate(corr =   corr, v = 0.1)
   
   
   lf[[l]] <- dae::rmvnorm(mean = rep(0,R),V = cov) #functional effect for effect l
 }
 
 
 # Plot the correlation matrix
 
 df_plot_block  <- data.frame(y= lf[[1]], x= 1:length(lf[[1]])) 
 
Pf_block <-  ggplot(df_plot_block, aes( x=x,y=y))+
   geom_line()+theme_void() +  # Minimal theme
   theme(   # Rotate x axis texts
     axis.title.x = element_blank(),
     axis.title.y = element_blank(),
     axis.text.x=element_blank(),
     axis.ticks.x=element_blank(),
     axis.text.y=element_blank(),
     axis.ticks.y=element_blank() ,
     legend.position = "none")
Pf_block
 library(ggplot2)
 library(reshape2)
 
 cor_matrix <- corr
 melted_cor_matrix <- melt(cor_matrix)
 
 # Plotting
P_block <-  ggplot(data = melted_cor_matrix, aes(Var1, Var2, fill = value)) +
   geom_tile(color = "white") +  # Use tiles, with white borders
   scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                        midpoint = 0, limit = c(-1,1), space = "Lab", 
                        name="Pearson\nCorrelation") +
   theme_void() +  # Minimal theme
   theme(   # Rotate x axis texts
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank() ,
         legend.position = "none") +
   coord_fixed()  # Ensure the tiles are square
P_block

set.seed(1)
f2 <-  simu_IBSS_per_level(7)$sim_func

f3 <-  simu_IBSS_per_level(7)$sim_func
f1 <-  simu_IBSS_per_level(7)$sim_func

idx_lst <-  susiF.alpha:::gen_wavelet_indx(7)


f1_list <- list()
wd1_list <- list()

xwd1 <-list()
ywd1 <-list()
library(wavethresh)
library(ashr)
for (s in 1:log2(length(f1))){
  tt <- rep( 0, length(f1) )
  
  
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

df1b$dummy <- rep( 0, nrow(df1b))
df1b <- df1b[-which( df1b$scale %in% c(1,2,3)),]
P_wac <-  ggplot( df1b, aes(x=x, y=y, colour=color))+
  #geom_point(size=2.5)+
  geom_hline(yintercept = 0)+
  facet_wrap(.~scale, ncol=1,strip.position = "left",scale="free")+
  scale_color_manual(values= c( "#377eb8","#e41a1c") )+
  geom_segment(data = df1b, aes(x =x, xend = x, 
                                y = dummy, yend = y), 
               
               colour = "#377eb8", 
               size = 1.5)+
  ggtitle("")+
  ylab("")+
  theme_bw()+  
  theme(
    legend.position = "none",
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank()
   # strip.placement = "outside"
    ) +
  xlab("") 
P_wac
df_plot_wac  <- data.frame(y= f1, x= 1:length(f1)) 

Pf_wac <-  ggplot(df_plot_wac , aes( x=x,y=y))+
  geom_line()+theme_void() +  # Minimal theme
  theme(   # Rotate x axis texts
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank() ,
    legend.position = "none")
Pf_wac
