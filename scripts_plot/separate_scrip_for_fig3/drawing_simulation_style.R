library(fsusieR)
library(cowplot)
library(ggplot2)

source(paste( path ,"/scripts_plot/plot_effect_benchmark.R", sep=""), echo=FALSE)
grid_plot <- ggdraw()+
  
  draw_plot(Pf_wac       ,
            x = 0.01 , y = .5, width = .33, height = .5)+
  draw_plot(P_wac         ,
            x = .0, y = .0, width = .33, height = .5)+
  
  draw_plot(Pf_block       ,
            x = .33 ,y = .5, width = .33, height = .5)+
  draw_plot(P_block         ,
            x = .33, y = .0, width = .33, height = .5)+
  draw_plot(Pf_decay       ,
            x = .66 , y = .5, width = .33, height = .5)+
  draw_plot(P_decay         ,
            x = .66, y = .0, width = .33, height = .5)
grid_plot




save_path=  paste0(getwd(),
                   "/plot/fig3_separate_panel/"
)
ggsave(grid_plot , file=paste0(save_path,"simulation_sketch.pdf"),
       width = 29.7,
       height = 21,
       units = "cm"
)