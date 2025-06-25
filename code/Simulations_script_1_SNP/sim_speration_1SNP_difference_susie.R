rm( list=ls())

# Load the required package
library(MASS)
library(fsusieR)
library(susieR)
lev_res= 7
effect_size=0.5
sd_noise=1


sim_function= function (lev_res= 7,
                        effect_size=0.5,
                        sd_noise=1,
                        n_diff=1){

  # Set parameters
  n <- nrow(N3finemapping$X)         # number of samples
  #number of individual having a different snp


  # Simulate
  X1 <-  N3finemapping$X[ , sample(size=1, 1:ncol(N3finemapping$X))]
  X2=X1

  val = unique(X1)

  id= sample( size = n_diff, 1:length(X1))

  for ( k in id){

    val_id=  which ( !(val==X2[k] ))
    X2[k]= val[ sample(size=1, val_id)]
  }

  # plot(X1,X2)
  X=cbind(X1,X2)
  # Optional: Check empirical correlation







  lf= rep(0,2^lev_res)
  lf[20:23]=effect_size

  dl= list()
  for ( i in 1:n){
    dl[[i]] = X[i,1]*lf + rnorm (n=2^lev_res, sd=sd_noise)
  }

  Y =do.call(rbind,dl)




    res= susie(y=Y[,21],X=X)

  tpv1= c()

    tt= summary(lm(Y[,21]~X[,1]))$coefficients[2,4]
    tpv1= c(tpv1,tt)



  tpv2= c()
    tt= summary(lm(Y[,21]~X[,2]))$coefficients[2,4]
    tpv2= c(tpv1,tt)



  pv1= min(tpv1)

  pv2= min(tpv2)


  out= list( cs=res$sets$cs$L1,
             pv1=pv1,
             pv2=pv2)
  return(out)
}

#out
#plot(-log10(tpv2), pch=19)
#points(-log10(tpv1), pch=19, col="green")


res=list()

teffect_size= seq(0.01,  1, length.out=100)
h=1
for ( o in teffect_size){
  for ( l in 1 :100 ){
    print(o)
    res[[h]]=sim_function(effect_size = o)

    h=h+1
  }


}






save(res, file="/project2/mstephens/fungen_xqtl/case_study/results/test_one_SNP_difference_susie.RData")


load("/project2/mstephens/fungen_xqtl/case_study/results/test_one_SNP_difference_susie.RData")
tt= list()
for (k in 1:length(res)){


  if(is.null(res[[k]]$cs)){
    tt[[k]]= c(0, res[[k]]$pv1)
  }else{
    if(length( res[[k]]$cs )==1 & ( res[[k]]$cs[[1]][1])==1){


      tt[[k]]= c(1, res[[k]]$pv1)
    }else{

      tt[[k]]= c(0, res[[k]]$pv1)
    }
  }



}

temp =do.call(rbind, tt)
plot( -log10(temp[,2] ), temp[,1])



df <- data.frame(response = temp[,1], predictor =-log10( temp[,2]))

# Step 2: Fit logistic regression model
model <- glm(response ~ predictor, data = df, family = binomial)

# Step 3: Create new data for smooth prediction curve
x_new <- data.frame(predictor = seq(min(df$predictor), max(df$predictor), length.out = 100))
x_new$predicted_prob <- predict(model, newdata = x_new, type = "response")

# Step 4: Plot with ggplot2
library(ggplot2)

P1= ggplot(df, aes(x = predictor, y = response)) +
  geom_point(size = 2) +
  geom_line(data = x_new, aes(x = predictor, y = predicted_prob), color = "blue", size = 1) +
  labs(title = "Logistic Regression Fit selection correct SNP only", x = "-log10 pv", y = "Probability seperating SNP") +
  theme_minimal()

P1


# Create breaks and labels
breaks <- seq(10, 60, by = 5)
labels <- paste0(breaks[-length(breaks)], "-", breaks[-1])

# Cut the predictor into bins
df$bin <- cut(df$predictor, breaks = breaks, right = FALSE, labels = labels)

# Compute mean response in each bin
bin_means <- aggregate(response ~ bin, data = df, FUN = mean, na.rm = TRUE)

print(bin_means)

library(ggplot2)
ggplot(bin_means, aes(x = bin, y = response)) +
  geom_col(fill = "steelblue") +
  labs(x = "-log10(pv) bin", y = "Mean response", title = "Mean trait per bin") +
  theme_minimal()



tt= list()
for (k in 1:length(res)){


  if(length( res[[k]]$cs[[1]])==1 & ( res[[k]]$cs[[1]][1])==2 ){


    tt[[k]]= c(1, res[[k]]$pv1)
  }else{

    tt[[k]]= c(0, res[[k]]$pv1)
  }

}

temp =do.call(rbind, tt)
plot( -log10(temp[,2] ), temp[,1])



df <- data.frame(response = temp[,1], predictor =-log10( temp[,2]))

# Step 2: Fit logistic regression model
model <- glm(response ~ predictor, data = df, family = binomial)

# Step 3: Create new data for smooth prediction curve
x_new <- data.frame(predictor = seq(min(df$predictor), max(df$predictor), length.out = 100))
x_new$predicted_prob <- predict(model, newdata = x_new, type = "response")

# Step 4: Plot with ggplot2
library(ggplot2)

P2=ggplot(df, aes(x = predictor, y = response)) +
  geom_point(size = 2) +
  geom_line(data = x_new, aes(x = predictor, y = predicted_prob), color = "blue", size = 1) +
  labs(title = "Logistic Regression Fit selection wrong SNP only", x = "-log10 pv", y = "Probability selecting the wrong  SNP only ") +
  theme_minimal()





library(gridExtra)
grid.arrange(P1,P2, ncol=2)

