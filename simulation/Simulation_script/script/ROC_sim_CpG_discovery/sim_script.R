



sim_perf_finding_CpG= function(n=100 ,
                               n_effect=5,
                               start_pos=10,
                               n_CPG=128,
                               h2=0.05){
  
  
  X= matrix(rnorm(n), ncol=1)
  
  noise_sd=   sqrt( ( 1/h2 )-1 )
  
  
  
  Y = matrix(rnorm(n*n_CPG,sd= noise_sd), ncol= n_CPG)
  
  Y[,start_pos:(start_pos+n_effect-1)]= Y[,start_pos:(start_pos+n_effect-1 )]+
    matrix(rep(X,n_effect), ncol =n_effect) 
  
  
  
  
  pv= do.call( c , lapply( 1:ncol(Y), function(j){
    summary(lm(Y[,j]~X))$coefficients[2,4]
  }))
  
  hmm_res= fsusieR:::univariate_HMM_regression(Y,X)
  
  TI_res= fsusieR:::univariate_TI_regression(Y,X)
  
  
  compute_pvalue <- function(estimate, lower_ci, upper_ci) {
    # Compute the standard error
    se <- (upper_ci - lower_ci) / (2 * qnorm(0.975))
    
    # Compute the test statistic
    z_score <- abs(estimate / se)
    
    # Compute the two-tailed p-value
    p_value <- 2 * (1 - pnorm(z_score))
    
    return(p_value)
  }
  
  # Example usage
  estimate <- 2.5
  lower_ci <- 1.0
  upper_ci <- 4.0
   
  
  
  pv_ti =do.call( c, lapply(  1: length(TI_res$effect_estimate),
                              function( i){
                                compute_pvalue(estimate=TI_res$effect_estimate[i], 
                                               lower_ci= TI_res$cred_band[2,i], 
                                               upper_ci= TI_res$cred_band[1,i])
                                
                              }
  )
  )
  
  pos_up <-  which(TI_res$cred_band [1,]<0)
  pos_low <- which(TI_res$cred_band [2,]>0)
  
  affected_TI= rep( 0,  length(TI_res$effect_estimate))
  if( length(pos_up)>0){
    affected_TI[pos_up]=1
  }
  if(length(pos_low)>0){
    affected_TI[pos_low]=1
  }
  CpG= rep(0, ncol(Y))
  CpG [start_pos:(start_pos+n_effect-1)]=1
  
  
  out = data.frame(CpG=CpG, pv=pv,
                   hmm_lfsr=hmm_res$lfsr,
                   pv_ti=pv_ti,
                   affected_TI=affected_TI,
                   n= rep( n, length(CpG)),
                   h2= rep( h2, length(CpG)))
  return( out)
  
}


#plot(-log10(pv_ti))
 