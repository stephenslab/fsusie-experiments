path <- getwd()
load(paste( path,"/simulation/Simulation_results/comparison_susie_fusie_128_sd1.RData", 
            sep=""))

 

tl =list()
for ( k in 1:length(res)){
  Number_effect = length(res[[k]]$true_pos)
  n_cs      = length(res[[k]]$susie_cs$cs)
  
  if( n_cs>0){
    n_false_effect=Reduce("+", lapply( 1:n_cs, function(l){
      ifelse( length(which(res[[k]]$true_pos%in%res[[k]]$susie_cs$cs[[l]] ))==0, 1,0)
    }))
  }else{
    n_false_effect=0 
  }
  n_effect= n_cs- n_false_effect
  tl[[k]]= c( Number_effect, n_cs, n_effect, n_false_effect)
}

df_susie= data.frame(do.call( rbind, tl))
colnames(df_susie)=c( "Number_effect", "n_cs", "n_effect", "n_false_effect")




