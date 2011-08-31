mode.density <- function(x,lbound=min(x,na.rm=TRUE),ubound=max(x,na.rm=TRUE), smooth=0.35, h=0){
#     posdensN <- locfit(~ x, alpha=c(2*smooth,0.3))
#     posdensN <- locfit(~ x)
#     xp <- floor(lbound):ceiling(ubound)
      x <- x[!is.na(x)]
      if(length(x)==0){return(NA)}
#     posdensN <- try(locfit( ~ lp(x, nn=2*smooth, h=h, maxk=500)),silent=TRUE)
#     if(inherits(posdensN,"try-error")){
       posdensN <- density(x, from=lbound, to=ubound)
#     }else{
#      xp <- seq(lbound, ubound, length=10000)
#      posdensN <- try(predict(posdensN, newdata=xp),silent=TRUE)
#      if(inherits(posdensN,"try-error")){
#       posdensN <- density(x, from=lbound, to=ubound)
#      }else{
#       posdensN <- list(x=xp,y=posdensN)
#      }
#     }
      posdensN$x[which.max(posdensN$y)]
}
ll.density <- function(x,lbound=min(x,na.rm=TRUE),ubound=max(x,na.rm=TRUE), smooth=0.35, h=0){
#     posdensN <- locfit(~ x, alpha=c(2*smooth,0.3))
#     posdensN <- locfit(~ x)
#     xp <- floor(lbound):ceiling(ubound)
      x <- x[!is.na(x)]
      if(length(x)==0){return(NA)}
      xp <- seq(lbound, ubound, length=10000)
      posdensN <- try(locfit( ~ lp(x, nn=2*smooth, h=h, maxk=500)),silent=TRUE)
      if(inherits(posdensN,"try-error")){
       posdensN <- density(x, from=lbound, to=ubound)
       list(x=xp,y=posdensN)
      }else{
       posdensN <- try(predict(posdensN, newdata=xp),silent=TRUE)
       if(inherits(posdensN,"try-error")){
        posdensN <- density(x, from=lbound, to=ubound)
        list(x=xp,y=posdensN)
       }else{
        list(x=xp,y=posdensN)
       }
      }
}
