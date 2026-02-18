load(paste0(getwd(),
            "/additional_work_for_revision/simulations/run_time.RData")
)
res

df=list()
for ( k in 1:length(res)){
  
df[[k]]=  c(res[[k]]$rt_fsusie[1],res[[k]]$rt_susie[1],2^ res[[k]]$s)
}
df=do.call(rbind,df)
df_plot= data.frame(run_time=  c(df[,1], df[,2]),
                    method= factor(rep(c("fSuSiE","SuSiE"),each = nrow(df))),
                    n_trait= factor(c(df[,3], df[,3])))
library(ggplot2)
ggplot(df_plot, aes( x=run_time, y=n_trait, colour = method))+
  geom_boxplot()






colors= c(  "darkblue","#D41159" ) #c(  "#FFC20A","#D41159" )
P1= ggplot(df_plot, aes( y=run_time, x=n_trait, colour = method))+
  geom_boxplot()+
  ylab("Run time")+
  xlab("Number of traits")+
  scale_color_manual(values=colors)+
  theme_linedraw()
P1
pdf(paste0 (path, "/plot/run_time_comp.pdf") , width = 10.27, height = 8.27)  
P1

dev.off()



load(paste0(getwd(),
            "/additional_work_for_revision/simulations/spike_run_time.RData")
)
res

df=list()
for ( k in 1:length(res)){
  
  df[[k]]=  c(res[[k]]$rt_fsusie[1],res[[k]]$rt_susie[1],2^ res[[k]]$s)
}
df=do.call(rbind,df)
df_plot= data.frame(run_time=  c(df[,1], df[,2]),
                    method= factor(rep(c("fSuSiE","SuSiE"),each = nrow(df))),
                    n_trait= factor(c(df[,3], df[,3])))
library(ggplot2)
ggplot(df_plot, aes( x=run_time, y=n_trait, colour = method))+
  geom_boxplot()










load(paste0(getwd(),
            "/additional_work_for_revision/simulations/run_time.RData")
)
res

df=list()
for ( k in 1:length(res)){
  
  df[[k]]=  c(res[[k]]$rt_fsusie[1],res[[k]]$rt_susie[1] ,2^ res[[k]]$s)
}
df=do.call(rbind,df)
df_plot= data.frame(run_time=  c(df[,1], df[,2]),
                    method= factor(rep(c("fSuSiE","SuSiE"),each = nrow(df))),
                    n_trait=  (c(df[,3], df[,3])))
library(ggplot2)
ggplot(df_plot, aes( y=run_time, x =n_trait, colour = method))+
  geom_point()+
  geom_smooth(method="lm")
lm(df[,2]~df[,3])
lm(df[,1]~df[,3])
