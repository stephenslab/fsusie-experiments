


n=100 
n_effect=5
start_pos=10
n_CPG=32
h2=0.05

X= matrix(rnorm(n), ncol=1)

noise_sd=   sqrt( ( 1/h2 )-1 )



Y = matrix(rnorm(n*n_CPG,sd= noise_sd), ncol= n_CPG)

Y[,start_pos:(start_pos+n_effect-1)]= Y[,start_pos:(start_pos+n_effect-1 )]+
                                      matrix(rep(X,n_effect), ncol =n_effect) 



 
pv= do.call( c , lapply( 1:ncol(Y), function(j){
  summary(lm(Y[,j]~X))$coefficients[2,4]
}))

plot( -log10(pv))
