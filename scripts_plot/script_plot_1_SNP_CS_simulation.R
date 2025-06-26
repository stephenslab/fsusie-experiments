 
library(susieR)
library(fsusieR)
library(mvsusieR)
library(gridExtra)
set.seed(123)
n_diff=1 # number of individual in which the SNP differs
# Set parameters
effect_size = 0.8
lev_res     = 7
sd_noise    = 1

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

plot(X1,X2)
X=cbind(X1,X2)
 

lf= rep(0,2^lev_res)
lf[20:30]=effect_size

dl= list()
for ( i in 1:n){
  dl[[i]] = X[i,1]*lf + rnorm (n=2^lev_res, sd=sd_noise)
}

Y =do.call(rbind,dl)
 
tpv1= c()
for ( j in 1:ncol(Y)){
  tt= summary(lm(Y[,j]~X[,1]))$coefficients[2,4]
  tpv1= c(tpv1,tt)
  
}
idx_susie=which.min(tpv1)

res0= susiF(Y=Y,X=X)

res1=  susie(y=Y[,idx_susie], X=X)

Y_t=Y[,20:30]
prior <- create_mixture_prior(R = ncol(  Y_t))
res2= mvsusie(X=X, Y= Y_t, prior_variance = prior)
 
res0$cs
res1$sets$cs
res2$sets$cs
 

 
load("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/outputs/1SNP_simulations/test_one_SNP_difference_susie_strongest_assoc.RData")
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

susie_res =do.call(rbind, tt)

df <- data.frame(response =susie_res[,1], predictor =-log10( susie_res[,2]))

# Step 2: Fit logistic regression model
model <- glm(response ~ predictor, data = df, family = binomial)

# Step 3: Create new data for smooth prediction curve
x_new <- data.frame(predictor = seq(min(df$predictor), max(df$predictor), length.out = 100))
x_new$predicted_prob <- predict(model, newdata = x_new, type = "response")

# Step 4: Plot with ggplot2
library(ggplot2)

Psusie= ggplot(df, aes(x = predictor, y = response)) +
  geom_point(size = 2) +
  geom_line(data = x_new, aes(x = predictor, y = predicted_prob), color = "blue", size = 1) +
  labs(title = "SuSiE", x = "-log10 pv", y = "Probability seperating SNP") +
  theme_minimal()


# Create breaks and labels
breaks <- seq(10, 60, by = 5)
labels <- paste0(breaks[-length(breaks)], "-", breaks[-1])

# Cut the predictor into bins
df$bin <- cut(df$predictor, breaks = breaks, right = FALSE, labels = labels)

# Compute mean response in each bin
bin_means <- aggregate(response ~ bin, data = df, FUN = mean, na.rm = TRUE)

print(bin_means)

library(ggplot2)
Psusiebin = ggplot(bin_means, aes(x = bin, y = response)) +
  geom_col(fill = "steelblue") +
  labs(x = "-log10(pv) bin", y = " Proportion of 1 SNP CS", title = "SuSiE")+
  theme_minimal()+ylim(c(0,1))


load("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/outputs/1SNP_simulations/test_one_SNP_difference_10cpgmvsusie.RData")
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

mvsusie_res=do.call(rbind, tt)



df <- data.frame(response = mvsusie_res[,1], predictor =-log10( mvsusie_res[,2]))

# Step 2: Fit logistic regression model
model <- glm(response ~ predictor, data = df, family = binomial)

# Step 3: Create new data for smooth prediction curve
x_new <- data.frame(predictor = seq(min(df$predictor), max(df$predictor), length.out = 100))
x_new$predicted_prob <- predict(model, newdata = x_new, type = "response")

# Step 4: Plot with ggplot2
library(ggplot2)

Pmvsusie = ggplot(df, aes(x = predictor, y = response)) +
  geom_point(size = 2) +
  geom_line(data = x_new, aes(x = predictor, y = predicted_prob), color = "blue", size = 1) +
  labs(title = "mvSuSiE", x = "-log10 pv", y = "Probability seperating SNP") +
  theme_minimal()


# Create breaks and labels
breaks <- seq(10, 60, by = 5)
labels <- paste0(breaks[-length(breaks)], "-", breaks[-1])

# Cut the predictor into bins
df$bin <- cut(df$predictor, breaks = breaks, right = FALSE, labels = labels)

# Compute mean response in each bin
bin_means <- aggregate(response ~ bin, data = df, FUN = mean, na.rm = TRUE)

print(bin_means)

library(ggplot2)
Pmvsusiebin = ggplot(bin_means, aes(x = bin, y = response)) +
  geom_col(fill = "steelblue") +
  labs(x = "-log10(pv) bin", y = " Proportion of 1 SNP CS", title = "mvSuSiE")+
  theme_minimal()+ylim(c(0,1))


load("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/outputs/1SNP_simulations/test_one_SNP_difference.RData")


tt= list()
for (k in 1:length(res)){
  
  
  if(length( res[[k]]$cs[[1]])==1 & ( res[[k]]$cs[[1]][1])==1){
    
    
    tt[[k]]= c(1, res[[k]]$pv1)
  }else{
    
    tt[[k]]= c(0, res[[k]]$pv1)
  }
  
}

temp =do.call(rbind, tt)

df <- data.frame(response = temp[,1], predictor =-log10( temp[,2]))

# Step 2: Fit logistic regression model
model <- glm(response ~ predictor, data = df, family = binomial)

# Step 3: Create new data for smooth prediction curve
x_new <- data.frame(predictor = seq(min(df$predictor), max(df$predictor), length.out = 100))
x_new$predicted_prob <- predict(model, newdata = x_new, type = "response")

# Step 4: Plot with ggplot2
library(ggplot2)

Pfsusie= ggplot(df, aes(x = predictor, y = response)) +
  geom_point(size = 2) +
  geom_line(data = x_new, aes(x = predictor, y = predicted_prob), color = "blue", size = 1) +
  labs(title = "fSuSiE", x = "-log10 pv", y = "Probability seperating SNP") +
  theme_minimal()


# Create breaks and labels
breaks <- seq(10, 60, by = 5)
labels <- paste0(breaks[-length(breaks)], "-", breaks[-1])

# Cut the predictor into bins
df$bin <- cut(df$predictor, breaks = breaks, right = FALSE, labels = labels)

# Compute mean response in each bin
bin_means <- aggregate(response ~ bin, data = df, FUN = mean, na.rm = TRUE)

print(bin_means)

library(ggplot2)
Pfsusiebin = ggplot(bin_means, aes(x = bin, y = response)) +
  geom_col(fill = "steelblue") +
  labs(x = "-log10(pv) bin", y = " Proportion of 1 SNP CS", title = "fSuSiE")+
  theme_minimal()+ylim(c(0,1))


Pout= grid.arrange(Psusie,Pmvsusie, Pfsusie, 
             ncol=3)
 

ggsave(Pout, file = "C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/plot/Logistic_performance.pdf", 
       width = 29.7,
       height = 21,
       units = "cm")

load("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/outputs/1SNP_simulations/test_one_SNP_difference_susie.RData")
tt= list()
for (k in 1:length(res)){
  
  
  if(is.null(res[[k]]$cs)){
    tt[[k]]= c(0, res[[k]]$pv1)
  }else{
    if(length( res[[k]]$cs )==1 & ( res[[k]]$cs[[1]][1])==2){
      
      
      tt[[k]]= c(1, res[[k]]$pv1)
    }else{
      
      tt[[k]]= c(0, res[[k]]$pv1)
    }
  }
  
  
  
}

susie_res =do.call(rbind, tt)

df <- data.frame(response =susie_res[,1], predictor =-log10( susie_res[,2]))

# Step 2: Fit logistic regression model
model <- glm(response ~ predictor, data = df, family = binomial)

# Step 3: Create new data for smooth prediction curve
x_new <- data.frame(predictor = seq(min(df$predictor), max(df$predictor), length.out = 100))
x_new$predicted_prob <- predict(model, newdata = x_new, type = "response")


nfalse= sum(susie_res[,1])
# Step 4: Plot with ggplot2
library(ggplot2)

Psusie= ggplot(df, aes(x = predictor, y = response)) +
  geom_point(size = 2) +
  geom_line(data = x_new, aes(x = predictor, y = predicted_prob), color = "blue", size = 1) +
  labs(title = paste0("SuSiE, false discovery = ",nfalse), x = "-log10 pv", y = "Probability selecting the wrong SNP") +
  theme_minimal()



load("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/outputs/1SNP_simulations/test_one_SNP_difference_10cpgmvsusie.RData")
tt= list()
for (k in 1:length(res)){
  
  
  if(is.null(res[[k]]$cs)){
    tt[[k]]= c(0, res[[k]]$pv1)
  }else{
    if(length( res[[k]]$cs )==1 & ( res[[k]]$cs[[1]][1])==2){
      
      
      tt[[k]]= c(1, res[[k]]$pv1)
    }else{
      
      tt[[k]]= c(0, res[[k]]$pv1)
    }
  }
  
  
  
}

mvsusie_res=do.call(rbind, tt)



df <- data.frame(response = mvsusie_res[,1], predictor =-log10( mvsusie_res[,2]))

# Step 2: Fit logistic regression model
model <- glm(response ~ predictor, data = df, family = binomial)

# Step 3: Create new data for smooth prediction curve
x_new <- data.frame(predictor = seq(min(df$predictor), max(df$predictor), length.out = 100))
x_new$predicted_prob <- predict(model, newdata = x_new, type = "response")

# Step 4: Plot with ggplot2
library(ggplot2)
nfalse= sum(mvsusie_res[,1])
Pmvsusie = ggplot(df, aes(x = predictor, y = response)) +
  geom_point(size = 2) +
  geom_line(data = x_new, aes(x = predictor, y = predicted_prob), color = "blue", size = 1) +
  labs(title =  paste0("mvSuSiE, false discovery = ",nfalse), x = "-log10 pv", y = "Probability selecting the wrong SNP") +
  theme_minimal()




load("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/outputs/1SNP_simulations/test_one_SNP_difference.RData")


tt= list()
for (k in 1:length(res)){
  
  
  if(length( res[[k]]$cs[[1]])==1 & ( res[[k]]$cs[[1]][1])==2){
    
    
    tt[[k]]= c(1, res[[k]]$pv1)
  }else{
    
    tt[[k]]= c(0, res[[k]]$pv1)
  }
  
}

temp =do.call(rbind, tt)

df <- data.frame(response = temp[,1], predictor =-log10( temp[,2]))

# Step 2: Fit logistic regression model
model <- glm(response ~ predictor, data = df, family = binomial)

# Step 3: Create new data for smooth prediction curve
x_new <- data.frame(predictor = seq(min(df$predictor), max(df$predictor), length.out = 100))
x_new$predicted_prob <- predict(model, newdata = x_new, type = "response")

# Step 4: Plot with ggplot2
library(ggplot2)
nfalse = sum(temp[,1])
Pfsusie= ggplot(df, aes(x = predictor, y = response)) +
  geom_point(size = 2) +
  geom_line(data = x_new, aes(x = predictor, y = predicted_prob), color = "blue", size = 1) +
  labs(title = paste0("fSuSiE, false discovery = ",nfalse), x = "-log10 pv", y = "Probability selecting the wrong SNP") +
  theme_minimal()

Pout2= grid.arrange(Psusie,Pmvsusie, Pfsusie, ncol=3)
 
ggsave(Pout2, file = "C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/plot/false_disc_performance.pdf", 
       width = 29.7,
       height = 21,
       units = "cm")
