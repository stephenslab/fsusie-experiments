rm(list=ls())

path <- getwd()
load(paste0(path,
            "/simulation/Simulation_script/script/additional_simualtion_for_fig3_panel_D_susie_pc_calibration_power/check_L_susie_128_distdecay_sd1.RData"))

 
 


df_susie= data.frame(do.call( rbind, res))
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


df= data.frame( T1= mean_T1_susie,
                Power= mean_power_susie,
                L= unique(df_simu$Number_effect))


plot( df$L, df$Power)
plot( df$L, df$T1)

