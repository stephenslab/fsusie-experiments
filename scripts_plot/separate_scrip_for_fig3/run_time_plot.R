

colors <- c( "#D41159","#1A85FF" )
##### run time ------

path <- getwd()
load(paste( path,"/simulation/Simulation_results/run_time_comp_p100.RData", 
            sep=""))

run_ncs_fsusie100    <- do.call( c, lapply( 1:length(res), function( i) length( res[[i]]$susiF_cs )))
run_ncs_sp_fsusie100 <- do.call( c, lapply( 1:length(res), function( i) length( res[[i]]$susiF_sp_cs)))

run_time_sp_fsusie100  <- do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_sp_time[3]))
run_time_fsusie100     <- do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_time[3]))


path <- getwd()
load(paste( path,"/simulation/Simulation_results/run_time_comp_p500.RData", 
            sep=""))
run_ncs_fsusie500    <- do.call( c, lapply( 1:length(res), function( i) length( res[[i]]$susiF_cs )))
run_ncs_sp_fsusie500 <- do.call( c, lapply( 1:length(res), function( i) length( res[[i]]$susiF_sp_cs)))


run_time_fsusie500     <- do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_time[3]))
run_time_sp_fsusie500  <- do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_sp_time[3]))

path <- getwd()
load(paste( path,"/simulation/Simulation_results/run_time_comp_p1000.RData", 
            sep="")) 

run_ncs_fsusie1000    <- do.call( c, lapply( 1:length(res), function( i) length( res[[i]]$susiF_cs )))
run_ncs_sp_fsusie1000 <- do.call( c, lapply( 1:length(res), function( i) length( res[[i]]$susiF_sp_cs)))

run_time_fsusie1000     <- do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_time[3]))
run_time_sp_fsusie1000  <- do.call( c, lapply( 1:length(res), function( i) res[[i]]$susiF_sp_time[3]))



df_run_time <- data.frame(runtime = c(run_time_fsusie100 ,
                                      run_time_fsusie500 ,
                                      run_time_fsusie1000,
                                      run_time_sp_fsusie100,
                                      run_time_sp_fsusie500,
                                      run_time_sp_fsusie1000),
                          ncs     = c(run_ncs_fsusie100,
                                      run_ncs_fsusie500,
                                      run_ncs_fsusie1000) ,
                          prior   = factor(c( rep("SPS",
                                                  (length(run_time_fsusie100)+length(run_time_fsusie500)+length(run_time_fsusie1000))),
                                              rep("IS",(length(run_time_sp_fsusie100)+length(run_time_sp_fsusie500)+length(run_time_sp_fsusie1000)))
                          )),
                          N= factor(c(rep( 100, (length(run_time_fsusie100) )),
                                      rep( 500, (length(run_time_fsusie500) )),
                                      rep( 1000,(length(run_time_fsusie1000) )),
                                      
                                      rep( 100, ( length(run_time_sp_fsusie100))),
                                      rep( 500, ( length(run_time_sp_fsusie500))),
                                      rep( 1000,( length(run_time_sp_fsusie1000)))
                                      
                          ))
                          
)
p_run_time <- ggplot(df_run_time, aes(runtime ,y= prior, col=prior))+
  geom_boxplot()+
  facet_wrap( N~. )+
  theme_linedraw() +
  ylab("Prior")+
  xlab("Run time (s)")+
  theme(strip.background =element_rect(fill="white"),
        strip.text = element_text(colour = 'black'))+
  scale_color_manual(values = colors)

p_run_time2 <- ggplot(df_run_time[ which( df_run_time$ncs <6 ),], aes(runtime ,y= as.factor(ncs), col=prior))+
  geom_boxplot()+
  facet_wrap(.~N  )+
  theme_linedraw() +
  ylab("Prior")+
  xlab("Run time (s)")+
  xlim(c(0,1000))+
  theme(strip.background =element_rect(fill="white"),
        strip.text = element_text(colour = 'black'),
        legend.position = "none")+
  scale_color_manual(values = colors)
p_run_time2

save_path=  paste0(getwd(),
                   "/plot/fig3_separate_panel/"
                   )
ggsave(p_run_time2 , file=paste0(save_path,"run_time.pdf"),
       width = 29.7,
       height = 10.5,
       units = "cm"
)
