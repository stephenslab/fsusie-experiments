library(susieRsmall)
library(fsusieR)
path_save = "D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/simulation/Simulation_script/script/additional_simualtion_for_fig3_panel_C_susie_pc_calibration_power"
res=list()
 
for (o  in (length(res)+1):10000) {
  
  L <- sample(1:20, size=1)#actual number of effect
  lf <-  list()
  for(l in 1:L){
    lf[[l]] <- simu_IBSS_per_level(lev_res=7)$sim_func #functional effect for effect l
  }
  
  list_files = list.files("D:/Document/Serieux/Travail/Package/1KG_data/1kg/rds")
  id = sample( 1:length(lf), size=1)
  G <- readRDS(paste0("D:/Document/Serieux/Travail/Package/1KG_data/1kg/rds/" ,list_files[id]))
  G <-   G[sample (1:nrow(  G), size=100, replace=FALSE), ]
  
  
   
  
  if( length(which(apply(G,2,var)==0))>0){
    G <- G[,-which(apply(G,2,var)==0)]
  }
  # G <- matrix( rnorm(100*300), nrow = 100)
  true_pos <- sample( 1:ncol(G), L)
  
  Y <- matrix(0 , ncol=  2^7 , nrow = 100)
  for ( i in 1:100){
    for ( l in 1:L){
      Y[i,] <- Y[i,]+ lf[[l]]*G[i,true_pos[[l]]]
    }
  }
  
  
  Y <-Y+matrix(rnorm((2^7)*100 ,sd=sd(c(Y))), nrow = 100)
  library(susieR)
  Z <- scale(Y,center = TRUE,scale = FALSE)
  out <- svd(Z)
  
  res <- susie(y=out$v[1,], X=G ,L=20)
  
   
  
  Number_effect = length( true_pos)
  n_cs      = length(res $cs)
  
  if( n_cs>0){
    n_false_effect=Reduce("+", lapply( 1:n_cs, function(l){
      ifelse( length(which( true_pos%in%res  $cs[[l]] ))==0, 1,0)
    }))
  }else{
    n_false_effect=0 
  }
  n_effect= n_cs- n_false_effect
  out <-   c( Number_effect, n_cs, n_effect, n_false_effect)
  
 
  res[[o]] <- (out)
  
  
  save(res, file=paste0(path_save, "/check_L_susie_128_sd1.RData"))
}
