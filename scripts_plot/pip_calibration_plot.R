path <- getwd()
load(paste( path,"/simulation/Simulation_results/comparison_susie_fusie_block_sd1.RData", 
            sep=""))


true_lab <- do.call( c,
                     lapply(1: length(res),
                            
                            function( i) {
                              
                              a <-  rep( 0,   length(res[[i]]$susiF_pip))
                              a[res[[i]]$true_pos] <- 1
                              return(a)
                            }
                            
                     )
)
data(N3finemapping)
X <- N3finemapping$X
True_cor <- cor(X)

 



score_fsusie <-  do.call( c, lapply( 1: length(res),
                                     function( i) res[[i]]$susiF_pip))
score_sp_fsusie <-  do.call( c, lapply( 1: length(res),
                                        function( i) res[[i]]$susiF_sp_pip))

score_susie <-  do.call( c, lapply( 1: length(res),
                                    function( i) res[[i]]$susie_rss_pip))


tdf= cbind(true_lab,score_sp_fsusie)
tdf=tdf[which(tdf[,1]==1),]
tdf 



set.seed(123)  # Set seed for reproducibility
n <- 10000       # Number of samples
data <- runif(n)  # Generate uniform random numbers

qqplot(qunif(ppoints(n)), sort(tdf[,2]),
       main = "QQ Plot - Uniform(0,1)",
       xlab = "Theoretical Quantiles",
       ylab = "Sample Quantiles",
       pch = 19, col = "blue")
abline(0, 1, col = "red", lwd = 2)  # Add 45-degree reference line


bin = c(0,0.1,0.2,0.3,0.4,.5,1)+0.0001
pip_bin= c(0)
y_axis= c()
 
errors =c(0)
for ( j in 1: (length(bin)-1)){
  
  n=length(tdf[,2])
  prop_est=length( which( tdf[,2] < bin [j+1]  ))/n
  pip_bin= c(pip_bin, prop_est
             )
  
  errors=  c(errors , sqrt((prop_est*(1-prop_est))/sqrt(n)))
  
  
}

 
plot( y=  pip_bin, x=bin, xlim=c(0,1),
      ylim=c(0,1))
points(x=bin,y=1.96*errors+  pip_bin, pch=20)
points(x=bin,y= -1.96*errors+  pip_bin, pch=20)

abline(a=0,b=1)






tdf= cbind(true_lab,score_susie)
tdf=tdf[which(tdf[,1]==1),]
tdf 



set.seed(123)  # Set seed for reproducibility
n <- 10000       # Number of samples
data <- runif(n)  # Generate uniform random numbers

qqplot(qunif(ppoints(n)), sort(tdf[,2]),
       main = "QQ Plot - Uniform(0,1)",
       xlab = "Theoretical Quantiles",
       ylab = "Sample Quantiles",
       pch = 19, col = "blue")
abline(0, 1, col = "red", lwd = 2)  # Add 45-degree reference line


bin = c(0,0.1,0.2,0.3,0.4,.5,1)+0.0001
pip_bin= c(0)
y_axis= c()

errors =c(0)
for ( j in 1: (length(bin)-1)){
  
  n=length(tdf[,2])
  prop_est=length( which( tdf[,2] < bin [j+1]  ))/n
  pip_bin= c(pip_bin, prop_est
  )
  
  errors=  c(errors , sqrt((prop_est*(1-prop_est))/sqrt(n)))
  
  
}


plot( y=  pip_bin, x=bin, xlim=c(0,1),
      ylim=c(0,1))
points(x=bin,y=1.96*errors+  pip_bin, pch=20)
points(x=bin,y= -1.96*errors+  pip_bin, pch=20)

abline(a=0,b=1)

