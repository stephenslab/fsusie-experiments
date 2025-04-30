rm(list=ls())

# Load necessary library
library(ggplot2)
library(dplyr)

save_path=  paste0(getwd(),
                   "/plot/fig3_separate_panel/"
)
path <- getwd()

## Gaussian ----
load(paste( path,"/simulation/Simulation_results/comparison_susie_fusie_128_sd1.RData", 
            sep=""))

 




size_fsusie <-  do.call( c, lapply( 1: length(res),
                                     function( i) lengths(res[[i]]$susiF_cs)))
size_sp_fsusie <-  do.call( c, lapply( 1: length(res),
                                        function( i) lengths(res[[i]]$susiF_sp_cs)))

size_susie <-  do.call( c, lapply( 1: length(res),
                                    function( i) lengths(res[[i]]$susie_cs$cs)))

idx1= which(size_fsusie> 2 & size_fsusie <6) 
idx3=which(size_fsusie> 10)
idx2=which(size_fsusie> 5 & size_fsusie <11)
size_fsusie[idx1] ="3-5"

size_fsusie[idx2] ="6-10"

size_fsusie[idx3] =">10"
table(size_fsusie)



idx1=which(size_sp_fsusie> 2 & size_sp_fsusie <6)
idx3=which(size_sp_fsusie> 10)
idx2= which(size_sp_fsusie> 5 & size_sp_fsusie <11)


c(table(size_sp_fsusie))

size_sp_fsusie[idx1] ="3-5"

c(table(size_sp_fsusie))
size_sp_fsusie[idx2] ="6-10"
c(table(size_sp_fsusie))
size_sp_fsusie[idx3] =">10"
c(table(size_sp_fsusie))


idx1=which(size_susie> 2 & size_susie <6)
idx3= which(size_susie> 10)
idx2= which(size_susie> 5 & size_susie <11)


size_susie[idx1] ="3-5"

size_susie[idx2] ="6-10"

size_susie[idx3] =">10"

c(table(size_fsusie))

c(table(size_sp_fsusie))
c(table(size_susie))


res_vec=c( c(table(size_fsusie)),
           
          c(table(size_sp_fsusie)),
          c(table(size_susie))
)

df_1= data.frame(
  res = res_vec,
  sim_type= rep( "Gaussian",length(res_vec) ),
  Method=  rep(c("fSuSiE SPS", "fSuSiE IS","SuSiE topPC"),each=5 ),
  cs_nb=   c(rep(c(">10","1","2","3-5", "6-10"),3)))
library(ggplot2)
library(dplyr)

# Make cs_nb a factor with correct order
df_1$cs_nb <- factor(df_1$cs_nb, levels = c("1", "2", "3-5", "6-10", ">10"))

# Plot
ggplot(df_1, aes(x = cs_nb, y = res, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge( )) +
  scale_fill_brewer(palette = "Set2", name = "Method") +
  xlab("Credible Set Size") +
  ylab("Number of Credible Sets") +
  ggtitle("Comparison of Credible Set Sizes across Methods") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12)
  )

##### Block ----
path <- getwd()
load(paste( path,"/simulation/Simulation_results/comparison_susie_fusie_block_sd1.RData", 
            sep=""))


size_fsusie <-  do.call( c, lapply( 1: length(res),
                                    function( i) lengths(res[[i]]$susiF_cs)))

c(table(size_fsusie))
size_sp_fsusie <-  do.call( c, lapply( 1: length(res),
                                       function( i) lengths(res[[i]]$susiF_sp_cs)))

size_susie <-  do.call( c, lapply( 1: length(res),
                                   function( i) lengths(res[[i]]$susie_cs$cs)))

idx1= which(size_fsusie> 2 & size_fsusie <6) 
idx3=which(size_fsusie> 10)
idx2=which(size_fsusie> 5 & size_fsusie <11)
size_fsusie[idx1] ="3-5"

size_fsusie[idx2] ="6-10"

size_fsusie[idx3] =">10"
table(size_fsusie)



idx1=which(size_sp_fsusie> 2 & size_sp_fsusie <6)
idx3=which(size_sp_fsusie> 10)
idx2= which(size_sp_fsusie> 5 & size_sp_fsusie <11)


c(table(size_sp_fsusie))

size_sp_fsusie[idx1] ="3-5"

c(table(size_sp_fsusie))
size_sp_fsusie[idx2] ="6-10"
c(table(size_sp_fsusie))
size_sp_fsusie[idx3] =">10"
c(table(size_sp_fsusie))


idx1=which(size_susie> 2 & size_susie <6)
idx3= which(size_susie> 10)
idx2= which(size_susie> 5 & size_susie <11)


size_susie[idx1] ="3-5"

size_susie[idx2] ="6-10"

size_susie[idx3] =">10"

c(table(size_fsusie))

c(table(size_sp_fsusie))
c(table(size_susie))


res_vec=c( 3.75*c(table(size_fsusie)),
           
           3.75*c(table(size_sp_fsusie)),
           3.75* c(table(size_susie))
)

df_2= data.frame(
  res = res_vec,
  sim_type= rep( "Block",length(res_vec) ),
  Method=  rep(c("fSuSiE SPS", "fSuSiE IS","SuSiE topPC"),each=5 ),
  cs_nb=   c(rep(c(">10","1","2","3-5", "6-10"),3)))

#### Decay ----


path <- getwd()
load(paste( path,"/simulation/Simulation_results/comparison_susie_fusie_distdecay_sd1.RData", 
            sep=""))

size_fsusie <-  do.call( c, lapply( 1: length(res),
                                    function( i) lengths(res[[i]]$susiF_cs)))
size_sp_fsusie <-  do.call( c, lapply( 1: length(res),
                                       function( i) lengths(res[[i]]$susiF_sp_cs)))

size_susie <-  do.call( c, lapply( 1: length(res),
                                   function( i) lengths(res[[i]]$susie_cs$cs)))


idx1= which(size_fsusie> 2 & size_fsusie <6) 
idx3=which(size_fsusie> 10)
idx2=which(size_fsusie> 5 & size_fsusie <11)
size_fsusie[idx1] ="3-5"

size_fsusie[idx2] ="6-10"

size_fsusie[idx3] =">10"
table(size_fsusie)



idx1=which(size_sp_fsusie> 2 & size_sp_fsusie <6)
idx3=which(size_sp_fsusie> 10)
idx2= which(size_sp_fsusie> 5 & size_sp_fsusie <11)


c(table(size_sp_fsusie))

size_sp_fsusie[idx1] ="3-5"

c(table(size_sp_fsusie))
size_sp_fsusie[idx2] ="6-10"
c(table(size_sp_fsusie))
size_sp_fsusie[idx3] =">10"
c(table(size_sp_fsusie))


idx1=which(size_susie> 2 & size_susie <6)
idx3= which(size_susie> 10)
idx2= which(size_susie> 5 & size_susie <11)


size_susie[idx1] ="3-5"

size_susie[idx2] ="6-10"

size_susie[idx3] =">10"

c(table(size_fsusie))

c(table(size_sp_fsusie))
c(table(size_susie))


res_vec=c( c(0,c(table(size_fsusie))),
           
           c(0,c(table(size_sp_fsusie))),
           c(table(size_susie))
)

df_3= data.frame(
  res = res_vec,
  sim_type= rep( "Decay",length(res_vec) ),
  Method= c(rep(c("fSuSiE SPS", "fSuSiE IS","SuSiE topPC"),each=5 )),
            
  cs_nb=   c(rep(c(">10","1","2","3-5", "6-10"),3)))



df= rbind(df_1,df_2,df_3)

df$ sim_type <- factor(df$ sim_type , levels = c("Gaussian", "Block", "Decay"))


P_out=ggplot(df , aes(x = cs_nb, y = res, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge( )) +
  scale_fill_brewer(palette = "Set2", name = "Method") +
  xlab("Credible Set Size") +
  ylab("Number of Credible Sets") +
  ggtitle("Comparison of Credible Set Sizes across Methods") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12)
  )+
  facet_grid(.~as.factor( sim_type))



save_path=  paste0(getwd(),
                   "/plot/cs_count"
)


ggsave(P_out , file=paste0(save_path,"cs_count.pdf"),
       width = 29.7,
       height = 21,
       units = "cm"
)
