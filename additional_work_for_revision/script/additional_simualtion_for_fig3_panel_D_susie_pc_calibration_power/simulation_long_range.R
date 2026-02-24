library(susieR)
library(dae)
library(fsusieR)
path_save = "C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/additional_work_for_revision/script/additional_simualtion_for_fig3_panel_D_susie_pc_calibration_power"
res=list()
library(fracdiff)
 
 


S=7

H=0.75
target_pve <- 0.20
for (o  in (length(res)+1):10000) {
  N=100 
  
  set.seed(o)
  list_files = list.files("C:/Document/Serieux/Travail/Package/1KG_data/1kg/rds")
  id = sample( 1:length(list_files), size=1)
  G <- readRDS(paste0("C:/Document/Serieux/Travail/Package/1KG_data/1kg/rds/" ,list_files[id]))
  G <-   G[sample (1:nrow(  G), size=100, replace=FALSE), ]
  
  
  if( length(which(apply(G,2,var)==0))>0){
    G <- G[,-which(apply(G,2,var)==0)]
  }
  
  L <- sample(1:10, size=1)#actual number of effect
  
  true_pos <- sample( 1:ncol(G), L)
  lf <-  list()
  for(l in 1:L){
    lf[[l]] <- simu_IBSS_per_level(lev_res=7)$sim_func #functional effect for effect l
  }
  R=2^7
  Y_signal <- matrix(0, nrow = 100, ncol = 2^S)
  
  for (i in 1:100) {
    for (l in 1:L) {
      Y_signal[i, ] <- Y_signal[i, ] + lf[[l]] * G[i, true_pos[l]]
    }
  }
  Y_noise <- matrix(0, nrow = 100, ncol = 2^S)
  
  for (i in 1:100) {
    eps <- fracdiff.sim(n = 2^S, d = H - 0.5)$series
    Y_noise[i, ] <- eps / sd(eps)
  }
  
  var_signal <- var(as.vector(Y_signal))
  var_noise  <- var(as.vector(Y_noise))
  
  
  signal_scale <- sqrt(
    (target_pve / (1 - target_pve)) * (var_noise / var_signal)
  )
  Y <- signal_scale * Y_signal + Y_noise
  
  var_signal <- var(as.vector(Y_signal))
  var_noise  <- var(as.vector(Y_noise))
  
  
  
  signal_scale <- sqrt(
    (target_pve / (1 - target_pve)) * (var_noise / var_signal)
  )
  Y <- signal_scale * Y_signal + Y_noise
  
  
  m1 <- susiF(Y=Y, X=G,L=20 ,L_start=11 ,nullweight= 10,  prior="mixture_normal", cal_obj =FALSE,  maxit=10)
  m2 <- susiF(Y=Y, X=G,L=20,L_start=11 ,nullweight= 10,post_processing = "none" ,
              prior="mixture_normal_per_scale" ,  maxit=10)
  
  PCA <- svd(Y)
  
  res_susie <-susieR::susie(X=G,
                            y=PCA$u[,1],
                            L=20
  ) 
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
  susie_res <-   c( Number_effect, n_cs, n_effect, n_false_effect,purity,cs_size)
  out= list( pip_fsusie_ps= m1$pip,
             pip_fsusie_nps= m2$pip,
             pip_susie= res_susie$pip,
             susie_res=susie_res,
             true_pos=true_pos
  )
  res[[o]] <- (out)
  
  save(res, file=paste0(path_save, "/long_128_.RData"))
}
