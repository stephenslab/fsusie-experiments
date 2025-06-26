 
path = paste0(getwd(),"/data/fig_4_data/")
res <- readRDS(paste0(path, "ROSMAP_mQTL.chr6_44880827_48309905.fsusie_mixture_normal_top_pc_weights.rds"))

fsusie_obj_me = res$`chr6:44880827-48309905`$ROSMAP_DLPFC_mQTL$fsusie_result
 
res_me <- readRDS(paste0(path, "ROSMAP_mQTL.chr6_44880827_48309905.fsusie_mixture_normal_top_pc_weights.input_data_1.rds"))
## to work from here



snp_names=attr( fsusie_obj_me$pip, "names")
pos_SNP_me <- as.numeric(sub("chr[0-9XY]+:([0-9]+):.*", "\\1", snp_names))
Y= as.data.frame(res_me$residual_Y)


X=as.data.frame(res_me$residual_X)
pos = as.data.frame(res_me$Y_coordinates) 
# this the causal SNP 
fsusie_obj_me$cs[[13]]


temp=list()
temp_pv=list()

temp_z=list()
for ( k in 1:ncol(Y)){
  
  fit_lm= lm(Y[,k]~X[,fsusie_obj_me$cs[[13]]])
   temp[[k]] = (summary(fit_lm)$coefficients[2,1] )*Y[,k]
   temp_pv[[k]] = (summary(fit_lm)$coefficients[2,4] ) 
   temp_z[[k]] = (summary(fit_lm)$coefficients[2,3] ) 
}


pv= do.call(c, temp_pv)
zscore= do.call(c, temp_z)
plot(-log10(pv))
plot(zscore)

#synthetic based on Beta
nY =Reduce("+",temp)
library(susieR)
m1= susie(y=nY, X=as.matrix(X))


temp=list()
for ( k in 1:ncol(Y)){
  
  fit_lm= lm(Y[,k]~X[,fsusie_obj_me$cs[[13]]])
  temp[[k]] = (summary(fit_lm)$coefficients[2,1]/summary(fit_lm)$coefficients[2,2] )*Y[,k]
  
}



#synthetic based on z score
nY2 =Reduce("+",temp)
library(susieR)
m2= susie(y=nY2, X=as.matrix(X))
m2$sets
# the fsusie result is 
fsusie_obj_me$cs[[13]]
cor( X[ c(m2$sets$cs[[1]] , fsusie_obj_me$cs[[13]])])

m1$sets
plot(X[,12183], X[, fsusie_obj_me$cs[[13]]])





