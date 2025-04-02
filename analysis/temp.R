bins <- c(0,1,2,5,10,20,Inf)
cs_size_susie  <- as.numeric(table(ha_snps_susie$cs))
cs_size_fsusie <- as.numeric(table(ha_snps_fsusie$cs))
cs_size_susie  <- cut(cs_size_susie,bins)
cs_size_fsusie <- cut(cs_size_fsusie,bins)
levels(cs_size_susie) <- bins[-1]
levels(cs_size_fsusie) <- bins[-1]
p1 <- ggplot(data.frame(cs_size = cs_size_susie),aes(x = cs_size)) +
  geom_histogram(stat = "count",color = "white",fill = "tomato",
                 width = 0.65) +
  labs(x = "CS size",y = "number of CSs",title = "SuSiE-topPC") +
  theme_cowplot(font_size = 10) +
  theme(plot.title = element_text(size = 10,face = "plain"))
p2 <- ggplot(data.frame(cs_size = cs_size_fsusie),aes(x = cs_size)) +
  geom_histogram(stat = "count",color = "white",fill = "tomato",
                 width = 0.65) +
  labs(x = "CS size",y = "number of CSs",title = "fSuSiE") +
  theme_cowplot(font_size = 10) +
  theme(plot.title = element_text(size = 10,face = "plain"))
plot_grid(p1,p2,nrow = 1,ncol = 2)
