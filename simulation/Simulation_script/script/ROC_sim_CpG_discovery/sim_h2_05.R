

path=getwd()

source(paste0( path,
               "/simulation/Simulation_script/script/ROC_sim_CpG_discovery/sim_script.R")
)




res= list()
for(  o in     1:100){
  
  res[[o]]= sim_perf_finding_CpG(h2=0.05,
                                 n =200,
                                 n_effect = 5)
  
  save(res, file =paste0(paste0(path, "/simulation/Simulation_script/script/ROC_sim_CpG_discovery/h2_05_n200_CpG_5.RData" )))
  print( o)
  
}






res= list()
for(  o in     1:100){
  
  res[[o]]= sim_perf_finding_CpG(h2=0.05,
                                 n =200,
                                 n_effect = 10)
  
  save(res, file =paste0(paste0(path, "/simulation/Simulation_script/script/ROC_sim_CpG_discovery/h2_05_n200_CpG_10.RData" )))
  print( o)
  
}

res= list()
for(  o in     1:100){
  
  res[[o]]= sim_perf_finding_CpG(h2=0.05,
                                 n =100,
                                 n_effect = 5)
  
  save(res, file =paste0(paste0(path, "/simulation/Simulation_script/script/ROC_sim_CpG_discovery/h2_05_n100_CpG_5.RData" )))
  print( o)
  
}



res= list()
for(  o in     1:100){
  
  res[[o]]= sim_perf_finding_CpG(h2=0.05,
                                 n =100,
                                 n_effect = 10)
  
  save(res, file =paste0(paste0(path, "/simulation/Simulation_script/script/ROC_sim_CpG_discovery/h2_05_n100_CpG_10.RData" )))
  print( o)
  
}

