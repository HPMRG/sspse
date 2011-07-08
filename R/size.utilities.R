mode.density <- function(x,lbound=min(x),ubound=max(x), smooth=0.35, h=0){
#     posdensN <- locfit(~ x, alpha=c(2*smooth,0.3))
#     posdensN <- locfit(~ x)
#     xp <- floor(lbound):ceiling(ubound)
      posdensN <- locfit( ~ lp(x, nn=2*smooth, h=h, maxk=500))
      xp <- seq(lbound, ubound, length=10000)
      posdensN <- list(x=xp,y=predict(posdensN, newdata=xp))
      posdensN$x[which.max(posdensN$y)]
}
ll.density <- function(x,lbound=min(x),ubound=max(x), smooth=0.35, h=0){
#     posdensN <- locfit(~ x, alpha=c(2*smooth,0.3))
#     posdensN <- locfit(~ x)
#     xp <- floor(lbound):ceiling(ubound)
      posdensN <- locfit( ~ lp(x, nn=2*smooth, h=h, maxk=500))
      xp <- seq(lbound, ubound, length=10000)
      list(x=xp,y=predict(posdensN, newdata=xp))
}
