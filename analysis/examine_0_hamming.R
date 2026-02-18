load("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/outputs/Hamming03.RData")


res= out$res

X=out$X
lt=list()
lt2=list()

id_1= which(lengths(res$cs)==1)
h=1
for ( l in id_1){
  
  x=X[,-res$cs[[l]]]
  tt= X[, (res$cs[[l]])]
  lt[[h]]=  max(cor(tt,x))
  lt2[[h]] =x[, which.max(cor(tt,x))]
  h=h+1
}
lt

h=6
l=id_1[h]
plot( X[, (res$cs[[l]])], lt2[[h]])

### CHR8 ----

load("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/outputs/Hamming03.RData")
res <- readRDS("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/outputs/ROSMAP_mQTL.chr8_63456198_66842322.fsusie_mixture_normal_top_pc_weights.rds")
input <- readRDS("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/outputs/ROSMAP_mQTL.chr8_63456198_66842322.fsusie_mixture_normal_top_pc_weights.rds.input_data.rds")

 
tt_c=  res$`chr8:63456198-66842322`$ROSMAP_DLPFC_mQTL$fsusie_result$cs
X=input$X0
lt=list()
lt2=list()

id_1= which(lengths(tt_c)==1)
h=1
for ( l in id_1){
  
  x=X[,-tt_c[[l]]]
  tt= X[, (tt_c[[l]])]
  lt[[h]]=  max(cor(tt,x))
  lt2[[h]] =x[, which.max(cor(tt,x))]
  h=h+1
}
lt


X=input$X0
  dim(X)
 
 length(res$`chr8:63456198-66842322`$ROSMAP_DLPFC_mQTL$fsusie_result$pip)
 
 
 
 res <- readRDS("../outputs/ROSMAP_mQTL.chr8_63456198_66842322.fsusie_mixture_normal_top_pc_weights.rds")
 input <- readRDS("../outputs/ROSMAP_mQTL.chr8_63456198_66842322.fsusie_mixture_normal_top_pc_weights.rds.input_data.rds")
 
 
 
 
 
 
 
 
 
 ##### CHR 11----- 
 

 res <- readRDS("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/outputs/ROSMAP_mQTL.chr11_77324757_82556425.fsusie_mixture_normal_top_pc_weights.rds")
 input <- readRDS("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/outputs/ROSMAP_mQTL.chr11_77324757_82556425.fsusie_mixture_normal_top_pc_weights.rds.input_data.rds")
 
 tt_c=  res$`chr11:77324757-82556425`$ROSMAP_DLPFC_mQTL$fsusie_result$cs
 X=input$X0
 lt=list()
 lt2=list()
 
 id_1= which(lengths(tt_c)==1)
 h=1
 for ( l in id_1){
   
   x=X[,-tt_c[[l]]]
   tt= X[, (tt_c[[l]])]
   lt[[h]]=  max(cor(tt,x))
   lt2[[h]] =x[, which.max(cor(tt,x))]
   h=h+1
 }
 lt
 
 
 X=input$X0
 
 dim(X)
 
 length(res$`chr11:77324757-82556425`$ROSMAP_DLPFC_mQTL$fsusie_result$pip)
 
 
 
 ##### CHR 6----- 
 
 
 res <- readRDS("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/outputs/ROSMAP_mQTL.chr6_7015726_11979320.fsusie_mixture_normal_top_pc_weights.rds")
 input <- readRDS("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/outputs/ROSMAP_mQTL.chr6_7015726_11979320.fsusie_mixture_normal_top_pc_weights.rds.input_data.rds") 
 
 tt_c=  res$`chr6:7015726-11979320`$ROSMAP_DLPFC_mQTL$fsusie_result$cs
 X=input$X0
 lt=list()
 lt2=list()
 
 id_1= which(lengths(tt_c)==1)
 h=1
 for ( l in id_1){
   
   x=X[,-tt_c[[l]]]
   tt= X[, (tt_c[[l]])]
   lt[[h]]=  max(cor(tt,x))
   lt2[[h]] =x[, which.max(cor(tt,x))]
   h=h+1
 }
 lt
 
 
 X=input$X0
 dim(X)
 
 length(res$`chr6:7015726-11979320`$ROSMAP_DLPFC_mQTL$fsusie_result$pip)
 
 
 