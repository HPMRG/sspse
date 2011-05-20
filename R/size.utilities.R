mode.density <- function(x,lbound=min(x),ubound=max(x), smooth=0.25){
      posdensN <- locfit(~ x, alpha=c(2*smooth,0.3))
#     posdensN <- locfit(~ x)
      xp <- floor(lbound):ceiling(ubound)
      posdensN <- list(x=xp,y=predict(posdensN, newdata=xp))
      posdensN$x[which.max(posdensN$y)]
}
