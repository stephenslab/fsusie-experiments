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

df_simu=df_susie
df_simu$power  <- df_simu$n_effect /df_simu$Number_effect 
df_simu$t1  <- 1- df_simu$n_false_effect /(df_simu$n_effect +df_simu$n_false_effect )
 



mean_power_susie <- rep( NA, max(df_simu$Number_effect) ) 
mean_T1_susie <- rep( NA,  max(df_simu$Number_effect) )
 
which_L <- rep( NA,  max(df_simu$Number_effect) )

h <- 1
for ( i in unique(df_simu$Number_effect))
{
  
  
  mean_power_susie[h] <- mean(df_simu$power [which(df_simu$Number_effect ==i )] )
    
  mean_T1_susie   [h] <- mean(df_simu$t1 [which(df_simu$Number_effect ==i  )] )
  
  
  which_L       [h] <- i
  
  
  h <- h+1
  
  
  
  
}

mean_T1_susie
mean_power_susie
