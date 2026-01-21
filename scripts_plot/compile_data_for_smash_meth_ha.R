library(smashr)
library(tools)

# M example.
set.seed(1)
dat <- readRDS(paste("ROSMAP_mQTL.chr20_53859688_57519449",
                     "fsusie_mixture_normal_top_pc_weights",
                     "input_data.rds",sep = "."))
i <- 1:256
y <- dat$Y[[1]]
y <- y[1,i]
pos <- dat$Y_coordinates[[1]]
pos <- pos[i,"start"]
res <- smash.gaus(y,v.est = FALSE)
plot(y,pch = 20,col = "black",xlab = "CpG",ylab = "beta value")
lines(res,type = "l",col = "royalblue")
methyl_dat <- data.frame(pos = pos,beta = y)

# HA example.
set.seed(1)
load("CR1_CR2_obj.RData")
names(obj_plot$count_df) <- c("CC","CT","TT","pos")
pos <- obj_plot$count_df$pos/1e6
i   <- which(pos > 207.45 & pos < 207.75)
dat <- obj_plot$count_df[i,]
y   <- dat$TT
pos <- dat$pos
rows1 <- which(y >= 5)
rows2 <- which(y < 5)
rows2 <- sample(rows2,3000)
rows  <- sort(c(rows1,rows2))
y     <- y[rows]
pos   <- pos[rows]
i     <- seq(1,2^13)
y     <- y[i]
pos   <- pos[i]
res <- smash.poiss(y,post.var = FALSE)
plot(pos/1e6,y,pch = 20,col = "black",xlab = "pos (Mb)",ylab = "read count")
lines(pos/1e6,res,type = "l",col = "royalblue")
h3k27ac_dat <- data.frame(pos = pos,y = y)
