library(susieR)
library(fsusieR)
path_save = "D:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/simulation/Simulation_script/script/additional_simualtion_for_fig3_panel_D_susie_pc_calibration_power"
res=list()
 
for (o  in (length(res)+1):10000) {
  
  
  N=100 
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
  # G <- matrix( rnorm(N*300), nrow = N)
  true_pos <- sample( 1:ncol(G), L)
  
  Y <- matrix(0 , ncol=  2^7 , nrow = N)
  for ( i in 1:N){
    for ( l in 1:L){
      Y[i,] <- Y[i,]+ lf[[l]]*G[i,true_pos[[l]]]
    }
  }
  
  
  Y <-Y+matrix(rnorm((2^7)*N ,sd=sd(c(Y))/sqrt(1)), nrow = N)
  
  PCA <- svd(Y)
  
  print("susie done")
  res_susie <-susieR::susie(X=G,
                    y=PCA$u[,1],
                    L=20
  )
  print(res_susie$sets)
  print("susie done")
  
  
  Number_effect = length( true_pos)
  n_cs      = length(res_susie$sets$cs)
  purity= mean ( res_susie$sets$purity[,1])
  cs_size= mean(lengths(res_susie$set$cs))
  
  if( n_cs>0){
    n_false_effect=Reduce("+", lapply( 1:n_cs, function(l){
      ifelse( length(which( true_pos%in%res_susie$sets$cs[[l]] ))==0, 1,0)
    }))
  }else{
    n_false_effect=0 
  }
  n_effect= n_cs- n_false_effect
  out <-   c( Number_effect, n_cs, n_effect, n_false_effect,purity,cs_size)
  
 
  res[[o]] <- (out)
  
  print( do.call(rbind,res))
  save(res, file=paste0(path_save, "/check_L_susie_128_sd1.RData"))
}
