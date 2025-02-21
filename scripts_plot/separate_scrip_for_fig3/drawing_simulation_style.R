library(fsusieR)
library(cowplot)
library(ggplot2)

source(paste( path ,"/scripts_plot/plot_effect_benchmark.R", sep=""), echo=FALSE)
grid_plot <- ggdraw()+
  
  draw_plot(Pf_wac       ,
            x = 0.0  , y = .6, width = .33, height = .4)+ 
  draw_plot(Pf_block       ,
            x = .35 ,y = .6, width = .28, height = .4)+
  
  draw_plot(Pf_decay       ,
            x = .68 , y = .6, width = .28, height = .4)+
  draw_plot(P_wac         ,
            x = .0, y = .0, width = .33, height = .6)+
  
 
  draw_plot(P_block         ,
            x = .33, y = .0, width = .33, height = .6)+
  draw_plot(P_decay         ,
            x = .66, y = .0, width = .33, height = .6)
grid_plot




save_path=  paste0(getwd(),
                   "/plot/fig3_separate_panel/"
)
ggsave(grid_plot , file=paste0(save_path,"simulation_sketch.pdf"),
       width = 25.7,
       height = 10.5,
       units = "cm"
)
