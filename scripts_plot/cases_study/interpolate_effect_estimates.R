# This is used to interpolate the fSuSiE effect estimates for the
# "zoom-in" plots.
interpolate_effect_estimates <- function (dat, pos) {
  n <- length(pos)
  out <- data.frame(pos    = pos,
                    effect = 0,
                    up     = 0,
                    low    = 0)
  for (i in 1:n) {
    up  <- which(dat$pos >= pos[i])
    low <- which(pos[i] >= dat$pos)
    du  <- abs(pos[i] - dat$pos[min(up)])
    di  <- abs(pos[i] - dat$pos[max(low)])
    dupdi <- du + di
    out[i,"effect"] <- (1 - du/dupdi) * dat[min(up),"effect"] +
                       (1 - di/dupdi) * dat[max(low),"effect"]
    out[i,"up"]     <- (1 - du/dupdi) * dat[min(up),"up"] +
                       (1 - di/dupdi) * dat[max(low),"up"]  
    out[i,"low"]    <- (1 - du/dupdi) * dat[min(up),"low"] +
                       (1 - di/dupdi) * dat[max(low),"low"] 
  }
  return(out)
}

