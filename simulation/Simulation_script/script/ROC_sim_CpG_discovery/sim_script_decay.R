 


sim_perf_finding_CpG= function(n=100 ,
                               n_effect=5,
                               start_pos=10,
                               n_CPG=128,
                               off_diag_corr = 0.99,
                               decay =0.9,
                               h2=0.05){
  
  
  X= matrix(rnorm(n), ncol=1)
  off_diag_corr = 0.99
  decay =0.9
  
  region_mat= matrix(off_diag_corr, ncol=  n_effect, nrow= n_effect)
  for ( i in 1: nrow(region_mat)){
    for ( j in  1:ncol(region_mat)){
      tt <- region_mat[i,j]
      region_mat[i,j] <- tt*decay ^( abs(i-j))
      region_mat[j,i] <- tt*decay ^( abs(i-j))
      
    }
  }
  
  effect=dae::rmvnorm(mean = rep(0, n_effect),V = region_mat)
  
  
  Y = matrix(0, ncol= n_CPG, nrow= nrow(X))
  for ( j in 1:n_effect){
    Y[, (start_pos+j-1) ]=  effect[j]*X
  }
  
  
  sd_Y_effect=sd(c(Y[,start_pos:(start_pos+n_effect-1)]))
  noise_sd= sqrt( (sd_Y_effect^2/h2)  -sd_Y_effect^2)
  
  Y <-Y+matrix(rnorm(prod(dim(Y)) ,
                     sd=noise_sd), nrow = nrow(Y))
  
  sigmoid <- function( x)
  {
    out <- 1/(1+exp(-x))
    return(out)
  }
  Y <- sigmoid(Y)
  
  
  pv= do.call( c , lapply( 1:ncol(Y), function(j){
    summary(lm(Y[,j]~X))$coefficients[2,4]
  }))
  
  hmm_res= fsusieR:::univariate_HMM_regression(Y,X)
  
  TI_res= fsusieR:::univariate_TI_regression(Y,X)
  
  
  grid= seq(from=0.00001, to =0.99, by =0.00001)
  cred_band = matrix( 0, ncol= length(TI_res$effect_estimate),
                      nrow=2)
  tl= list()
  for( i in 1:length(grid)){
    alpha=grid[i]
    coeff= qnorm(1- alpha /2)
    
    cred_band=0*cred_band
    cred_band [1, ]= TI_res$effect_estimate +  coeff* sqrt(TI_res$fitted_var )
    cred_band [2, ]= TI_res$effect_estimate -  coeff* sqrt(TI_res$fitted_var )
    
    tout=rep( 0, ncol(cred_band))
    
    idxup = which( cred_band[1,]<0)
    idxlow= which( cred_band[2,]>0)
    if( length(idxup)>0){
      tout[idxup]=1
    }
    if( length(idxlow)>0){
      tout[idxlow]=1
    }
    tl[[i]]=tout
  }
  
  tmat= do.call(rbind, tl)
  
  thresh=   do.call(c, 
                    lapply( 1: ncol( tmat),
                            function( i ){
                              if(length( which ( tmat[,i]>0))==0){#ie only zero
                                idx= grid[nrow(tmat)]
                              }else{
                                idx= grid[ min( which ( tmat[,i]>0))]#highest  level of condfidence 
                                # in which the credible band do not include 0 
                              }
                              return(idx)
                            }
                    )
  )
  
  
  compute_pvalue <- function(estimate, lower_ci, upper_ci) {
    # Compute the standard error
    se <- (upper_ci - lower_ci) / (2 * qnorm(0.975))
    
    # Compute the test statistic
    z_score <- abs(estimate / se)
    
    # Compute the two-tailed p-value
    p_value <- 2 * (1 - pnorm(z_score))
    
    return(p_value)
  }
  
  pv_ti =do.call( c, lapply(  1: length(TI_res$effect_estimate),
                              function( i){
                                compute_pvalue(estimate=TI_res$effect_estimate[i], 
                                               lower_ci= TI_res$cred_band[2,i], 
                                               upper_ci= TI_res$cred_band[1,i])
                                
                              }
  )
  )
  
  #pos_up <-  which(TI_res$cred_band [1,]<0)
  #pos_low <- which(TI_res$cred_band [2,]>0)
  
  #affected_TI= rep( 0,  length(TI_res$effect_estimate))
  #if( length(pos_up)>0){
  #  affected_TI[pos_up]=1
  #}
  #if(length(pos_low)>0){
  #  affected_TI[pos_low]=1
  #}
  CpG= rep(0, ncol(Y))
  CpG [start_pos:(start_pos+n_effect-1)]=1
  
  
  out = data.frame(CpG=CpG, pv=pv,
                   hmm_lfsr=hmm_res$lfsr,
                   pv_ti=pv_ti,
                   affected_TI=thresh,
                   n= rep( n, length(CpG)),
                   h2= rep( h2, length(CpG)))
  return( out)
  
}


#plot(-log10(pv_ti))
