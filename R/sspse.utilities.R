#' @keywords internal
mode.density <- function(x,lbound=min(x,na.rm=TRUE),ubound=max(x,na.rm=TRUE), smooth=0.35, h=0){
#     posdensN <- locfit(~ x, alpha=c(2*smooth,0.3))
#     posdensN <- locfit(~ x)
#     xp <- floor(lbound):ceiling(ubound)
      x <- x[!is.na(x)]
      if(length(x)==0){return(NA)}
      if(length(x[x>=lbound&x<=ubound])==0){return(NA)}
      if(length(x[x>=lbound&x<=ubound])==1){return(x[x>=lbound&x<=ubound])}
#     posdensN <- try(locfit( ~ lp(x, nn=2*smooth, h=h, maxk=500)),silent=TRUE)
#     if(inherits(posdensN,"try-error")){
       posdensN <- stats::density(x, from=lbound, to=ubound)
#     }else{
#      xp <- seq(lbound, ubound, length=10000)
#      posdensN <- try(predict(posdensN, newdata=xp),silent=TRUE)
#      if(inherits(posdensN,"try-error")){
#       posdensN <- stats::density(x, from=lbound, to=ubound)
#      }else{
#       posdensN <- list(x=xp,y=posdensN)
#      }
#     }
      posdensN$x[which.max(posdensN$y)]
}
#' @keywords internal
ll.density <- function(x,lbound=min(x,na.rm=TRUE),ubound=max(x,na.rm=TRUE), smooth=0.35, h=0, method="bgk"){
#     posdensN <- locfit(~ x, alpha=c(2*smooth,0.3))
#     posdensN <- locfit(~ x)
#     xp <- floor(lbound):ceiling(ubound)
      x <- x[!is.na(x)]
      if(length(x)==0){return(NA)}
      xp <- seq(lbound, ubound, length=10000)
      if(method == "bgk"){
#       posdensN <- try(locfit( ~ lp(x, nn=2*smooth, h=h, maxk=500)),silent=TRUE)
        posdensN <- .catchToList(bgk_kde(x,n=2^(ceiling(log(ubound-lbound)/log(2))),MIN=lbound,MAX=ubound))
        if(!is.null(posdensN$error)){
         posdensN <- stats::density(x, from=lbound, to=ubound)
         list(x=xp,y=posdensN)
        }else{
#        posdensN <- try(predict(posdensN, newdata=xp),silent=TRUE)
         posdensN <- .catchToList(stats::spline(x=posdensN$value[1,],y=posdensN$value[2,],xout=xp)$y)
         if(!is.null(posdensN$error)){
          posdensN <- stats::density(x, from=lbound, to=ubound)
          list(x=xp,y=posdensN)
         }else{
          list(x=xp,y=posdensN$value)
         }
        }
      }else{
        if(requireNamespace("densEstBayes", quietly = TRUE)){
          control<-list(samples=4000, burnin=1000)
          a=densEstBayes::densEstBayes(x,method="NUTS",
            control=densEstBayes::densEstBayes.control(range.x=c(lbound,ubound),numBins=401,numBasis=round(smooth*50/0.35),
                                                       nKept=control$samples,nWarm=control$burnin))
          xTrang <- seq(-0.05, 1.05, length = length(xp))
          Xg <- cbind(1,xTrang)
          Zg <- .ZOSull(xTrang,intKnots=a$intKnots,range.x=c(-0.05,1.05))
          Cg <- cbind(Xg,Zg)
          betauMCMC <- a$stochaFitObj$betauMCMC
          etaHatMCMC <- crossprod(t(Cg),betauMCMC)
          posdensN <- exp(apply(etaHatMCMC, 1, mean))
          posdensN <- length(xp)*posdensN / ((lbound-ubound)*sum(posdensN))
          list(x=xp,y=posdensN)
        }
      }
}
#' @keywords internal
HPD.density <- function(x,lbound=min(x,na.rm=TRUE),ubound=max(x,na.rm=TRUE)){
      x <- x[!is.na(x)]
      if(length(x)==0){return(c(NA,NA))}
      posdensN <- stats::density(x, from=lbound, to=ubound)
      cy <- cumsum(posdensN$y/sum(posdensN$y))
      cy <- c(posdensN$x[which.max(cy>0.025)],
              posdensN$x[which.max(cy>0.975)])
      if(is.na(cy[1])) cy[1] <- lbound
      if(is.na(cy[2])) cy[2] <- ubound
      cy
}

.catchToList <- function(expr) {
  val <- NULL
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, w$message)
    invokeRestart("muffleWarning")
  }
  myError <- NULL
  eHandler <- function(e) {
    myError <<- e$message
    NULL
  }
  val <- tryCatch(withCallingHandlers(expr, warning = wHandler), error = eHandler)
  list(value = val, warnings = myWarnings, error=myError)
} 

# Compute the horvitz thompson estimate
HT.estimate<-function(weights,outcome){
  nas <- is.na(weights)|is.na(outcome)
  if(is.factor(outcome)){
    num <- rep(0,length=length(levels(outcome[!nas])))
    a <- tapply(weights[!nas],as.numeric(outcome[!nas]),sum)
    num[as.numeric(names(a))] <- a
  }else{
    num<-sum(outcome[!nas]*weights[!nas])
  }
  den<-sum(weights[!nas])
  num/den
}
