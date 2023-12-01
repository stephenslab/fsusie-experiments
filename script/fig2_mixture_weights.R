library(ggplot2)
library(cowplot)
dat <- readRDS("../data/fig_2_data/fitted_weight_data.rds")
pdat <- data.frame(scale = rep(0:7,times = 32),
                   k     = rep(1:32,each = 8),
                   x     = as.vector(dat))
pdat <- transform(pdat,x = cut(x,c(-1,0.01,0.9,1)))
p <- ggplot(pdat,aes(x = k,y = scale,fill = x)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("white","dodgerblue","darkblue")) +
  scale_y_continuous(breaks = 0:7) +
  theme_cowplot(font_size = 9)
ggsave("pi_matrix.eps",p,height = 1.2,width = 4)
