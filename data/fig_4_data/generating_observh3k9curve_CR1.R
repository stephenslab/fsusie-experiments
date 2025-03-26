rm(list=ls(
  
))
library(data.table)

h3k9 <- readRDS("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/haQTL_CR1_raw_peaks.rds")

meta_data= fread("C:/Document/Serieux/Travail/Data_analysis_and_papers/xQTL_data/H3K9ac_count/H3K9acDomains.csv" )


Id_map = fread("C:/Document/Serieux/Travail/Data_analysis_and_papers/xQTL_data/H3K9ac_count/sampleSheetAfterQc.csv")


 
idcount= as.numeric(unique(h3k9$samples))
h3k9$IDsample=NA

l_name=rep(NA, length(unique(h3k9$samples)))
for (i in 1: length(idcount)){
  
   
  
  if ( length(Id_map$SampleID[which(Id_map$ProjID ==idcount[i])] )>0){
    h3k9$IDsample[ which( h3k9$samples==unique(h3k9$samples)[i])] =Id_map$SampleID[which(Id_map$ProjID ==idcount[i])]
  }
  
}
 


h3k9 = h3k9[-which(is.na(h3k9$IDsample)),]
h3k9$start <- as.numeric(h3k9$start)
h3k9$end <- as.numeric(h3k9$end)
h3k9$normalized_reading <- as.numeric(h3k9$normalized_reading)










data= h3k9


h3k9 [which(h3k9$start==h3k9$start[1]),]




plot(data$start[1:20000],
     data$normalized_reading [1:20000])






res_ha <- readRDS("c:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/fsusie_object/raw_data/ROSMAP_haQTL.chr1_205117782_208795513.fsusie_mixture_normal_top_pc_weights.input_data (1).rds")


X=as.data.frame(res_ha$residual_X) 


lead_SNP= matrix(X[,10765 ], ncol=1)

hist(lead_SNP)

binarized_SNP=0* lead_SNP 

for ( i in 1:nrow(lead_SNP)){
  if(lead_SNP[i]<0){
    binarized_SNP[i]=1
    if(lead_SNP[i]< -0.8){
       binarized_SNP[i]=2
    }
  }
}

table(binarized_SNP)


rownames(binarized_SNP)= rownames(X)

h3k9$SNP = rep(0,nrow(h3k9))

h3k9$SNP [which (h3k9$IDsample %in% row.names(binarized_SNP)[which(binarized_SNP==1)])]=1


h3k9$SNP [which (h3k9$IDsample %in% row.names(binarized_SNP)[which(binarized_SNP==2)])]=2
table(h3k9$SNP)

unique(h3k9$samples)


h3k9 [c(which(h3k9$samples =="00482428")[1:20],
        which(h3k9$samples =="00668310")[1:20]),]
        





plot(h3k9$normalized_reading [ which(h3k9$samples =="00668310")][1:2000], type="l")
plot(h3k9$normalized_reading [ which(h3k9$samples =="00482428")][1:2000], type="l")








start = min(h3k9$start)


end= max(h3k9$end)

obs_pos= start:end


mean_func2= 0*obs_pos
idx = unique( h3k9$IDsample[which (      h3k9$SNP==2)])
for ( i in 1: length(idx )){
  
  sub= h3k9[which(h3k9$IDsample==idx  [i]),]
  
  for ( k in 1:nrow(sub)){
    
    sub_pos = which( obs_pos>=  sub$start[k ] & obs_pos  <= sub$end[k ])
    if (length(sub_pos)>0){
      
      mean_func2[ sub_pos]= mean_func2[ sub_pos]+sub$normalized_reading[k]
    }
    
  }
  
  print(i)
  
} 


mean_func2=mean_func2/length(idx )






mean_func1= 0*obs_pos
idx = unique( h3k9$IDsample[which (      h3k9$SNP==1)])
for ( i in 1: length(idx )){
  
  sub= h3k9[which(h3k9$IDsample==idx  [i]),]
  
  for ( k in 1:nrow(sub)){
    
    sub_pos = which( obs_pos>=  sub$start[k ] & obs_pos  <= sub$end[k ])
    if (length(sub_pos)>0){
      
      mean_func1[ sub_pos]= mean_func1[ sub_pos]+sub$normalized_reading[k]
    }
    
  }
  
  print(i)
  
}
mean_func1=mean_func1/length(idx )



mean_func0= 0*obs_pos
idx = unique( h3k9$IDsample[which (      h3k9$SNP==0)])
for ( i in 1: length(idx )){
  
  sub= h3k9[which(h3k9$IDsample==idx  [i]),]
  
  for ( k in 1:nrow(sub)){
    
    sub_pos = which( obs_pos>=  sub$start[k ] & obs_pos  <= sub$end[k ])
    if (length(sub_pos)>0){
      
      mean_func0[ sub_pos]= mean_func0[ sub_pos]+sub$normalized_reading[k]
    }
    
  }
  
  print(i)
  
}
mean_func0=mean_func0/length(idx )



obs_curve= list(mean_func0=mean_func0,
                mean_func1=mean_func1,
                mean_func2=mean_func2,
                obs_pos=obs_pos)

save(obs_curve,
     file="C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/data/fig_4_data/row_count_CR1.RData")
 