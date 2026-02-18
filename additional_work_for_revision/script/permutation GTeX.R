 
set.seed(1)
 path= paste0(getwd(),"/additional_work_for_revision/data_analysis/GTeX/" )
lf=list.files(path)
tt1= list()
tt2=list()
tt1= list()
tt0=list()
for (k in (length(tt1)+1):100){
  id = sample (1: length(lf))[1]
  load(paste0(path , lf[id]))
  X_perm=as.matrix(   out$X[ sample(1:nrow(out$X)),])
  tt= fsusieR::susiF(Y=log1p(out$Y), X=X_perm, L=2, nullweight =.1, post_processing = "none" , quantile_trans = TRUE,
                     thresh_lowcount=.4
  )
  lol=susieR::susie(y=log1p(apply(out$Y/out$size_factor_local,1, sum)), X=X_perm)
  tt1[[k]]=tt$cs
  tt= fsusieR::susiF(Y=log1p(out$Y), X=X_perm, L=2, nullweight =.1, post_processing = "none" , quantile_trans = TRUE,cor_small = TRUE,
                     thresh_lowcount=.4 
  )
  tt0[[k]]=tt$cs
  tt2[[k]]=lol$sets
  print(
    
    lapply( 1:length(tt1), function(i)length(tt1[[i]][[1]] )))
  print(
    
    lapply( 1:length(tt0), function(i)length(tt0[[i]][[1]] )))
}


out= list(tt1=tt1,
          tt2=tt2)
save()
 plot(tt$pip)
lengths(tt1)

print(
  lapply( 1:length(tt1), function(i)length(tt2[[i]]$cs)))

lapply( 1:length(tt1), function(i)length(tt1[[i]][[1]] ))
plot(tt$lBF[[1]])

plot(tt$fitted_func[[1]])awdawd