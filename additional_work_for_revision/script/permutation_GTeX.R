
set.seed(1)
path= paste0("/project2/mstephens/wdenault/GTEX_analysis_Fsusie/result_Gao/" )
lf=list.files(path)
tt1= list()
tt2=list()
tt0=list()
tt_param=list()

o =1
for ( null_w in c( 0.1, 0.5))
{
  for ( thres_c in  c( 0.05,0.1, 0.5)){


    for (k in 1:10){
      id = sample (1: length(lf))[1]
      load(paste0(path , lf[id]))
      X_perm=as.matrix(   out$X[ sample(1:nrow(out$X)),])
      tt= fsusieR::susiF(Y=log1p(out$Y), X=X_perm, L=2,
                         nullweight =null_w,
                         post_processing = "none" ,
                         quantile_trans = TRUE,
                         thresh_lowcount=thres_c
      )


      tt_small= fsusieR::susiF(Y=log1p(out$Y), X=X_perm, L=2,
                         nullweight =null_w,
                         post_processing = "none" ,
                         quantile_trans = TRUE,
                         cor_small = TRUE,
                         thresh_lowcount=thres_c
      )
      lol=susieR::susie(y=log1p(apply(out$Y/out$size_factor_local,1, sum)), X=X_perm)
      tt0[[o]]=tt_small$cs
      tt1[[o]]=tt$cs

      tt2[[o]]=lol$sets

      tt_param[[o]]= c(null_w, thres_c)
      out= list(fsusie=tt1,
                fsusie_corsmall=tt0,
                susie=tt2,
                param= tt_param)

      save(out , file="/home/wdenault/fsusi_simu/correlated_sim/permutation_GTEX.RData")
      #print(

       # lapply( 1:length(tt1), function(i)length(tt1[[i]][[1]] )))

      #print(
       # lapply( 1:length(tt1), function(i)length(tt2[[i]]$cs)))
      o= o+1
    }


  }


}






