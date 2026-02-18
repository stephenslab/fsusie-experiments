library(bfg)


load("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/additional_work_for_revision/data_analysis/GTeX/APP.RData")


n = nrow(out$Y)
m =ncol(out$Y)
p =ncol(out$X) #200
p0 = 5
RSNR = 5
ell0 = 0.1
rho = 0.75

data_generated = gen_data(n,m,p,p0,RSNR,ell0,rho,re = T)

t = data_generated$T
tau0_prime0 = data_generated$tau0_prime
fit = bfg(out$Y,out$X,t,tau0_prime0,data_generated=data_generated,
          interactions = F, plotting = F,thinning = 1, N.iter=2000)
