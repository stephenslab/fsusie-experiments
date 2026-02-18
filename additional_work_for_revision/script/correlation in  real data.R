
path= getwd()
dat <- readRDS(paste0(path,
                      "/additional_work_for_revision/data_analysis/ROSMAP_mQTL.chr12_31679199_37240000.fsusie_mixture_normal_top_pc_weights.with_input.rds")
)

image(cor(dat$original_data$residual_Y$ROSMAP_DLPFC_mQTL))
image(cor(dat$original_data$Y[[1]]
          ))
Y= dat$original_data$residual_Y$ROSMAP_DLPFC_mQTL
map_data <- fsusieR:::remap_data(Y=Y,
                       pos=dat$original_data$Y_coordinates[[1]]$start,
                       max_scale=10)

outing_grid <- map_data$outing_grid
Y           <- map_data$Y
library(fsusieR)

image(cor(Y)) 
mean(cor(Y))
image(cor( Y- dat$fsusie_result$ind_fitted_func))

mean(cor( Y- dat$fsusie_result$ind_fitted_func))
X= matrix(dat$original_data$X_data[[1]], nrow=636)
res= susiF(Y=Y,X=X, L=3)



plot(res$fitted_func[[3]])
lines(res$cred_band[[3]][1,])
lines(res$cred_band[[3]][2,])
predicted= 0*Y
for ( k in 1 : length(res$cs)){
  
  predicted=predicted+ X[, res$cs[[k]][1]]%*%t(res$fitted_func[[k]])
  
}

image(cor( Y-predicted  ))
image(cor( Y   ))


mean(cor(Y))

mean(cor(Y-predicted  ))

plot(c(cor(Y)) -c(cor( Y-predicted  )))



res$cs


pv= c()

for ( k in 1:ncol(Y)){
  
  pv= c(pv,summary(lm(Y[,k]~X[,296]))$coefficients[2,4])
}
plot(-log10(pv))
