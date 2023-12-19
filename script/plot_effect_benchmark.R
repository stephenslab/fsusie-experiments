set.seed(4)

L <- sample(1:5, size =1)#actual number of effect
lf <-  list()
R=2^7
for(l in 1:L){
  
  block_cov_generate = function(corr=0.9, v=1){
    R = dim(corr)[1]
    cov = diag(sqrt(v),R) %*% corr %*% diag(sqrt(v),R)
    return(cov)
  }
  
  decay_corr_generate = function(R =2^7,# # sites
                                 num_block = 3,# #islands within the region
                                 off_diag_corr = 0.99,
                                 decay =0.99#
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
   theme_minimal() +  # Minimal theme
 
  theme_minimal() +  # Minimal theme
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
 
 plot( lf[[1]])
 image (cov)
 # Plot the correlation matrix
 plot( lf[[1]])
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
   theme_minimal() +  # Minimal theme
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



dat <- readRDS("B_mat.rds")
clrs <- c("white","#fee5d9","#fcae91","#fb6a4a","#de2d26","#a50f15")
image(abs(as.matrix(dat[40:80,100:128])),col = clrs)
 