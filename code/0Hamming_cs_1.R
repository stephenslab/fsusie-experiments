library(fsusieR)
dta <- readRDS("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/outputs/ROSMAP_mQTL.chr6_7015726_11979320.fsusie_mixture_normal_top_pc_weights.rds.input_data.rds")

Y= dta$Y1
X=dta$X1
res = susiF(Y=Y,X=X, L=20)
out=  list(data=dta,
           Y=Y,
           X=X,
           res=res)
save(out, file= "C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/outputs/Hamming01.RData")