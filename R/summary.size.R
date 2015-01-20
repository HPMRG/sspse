summary.psess <- function(object, support=1000, HPD.level=0.95, ...){
#summary.psess <- function(object, ...){
  p.args <- as.list( sys.call() )[-c(1,2)]
  formal.args<-formals(sys.function())[-1]
# control <- list(support=1000,HPD.level=0.95)

  control<-list()
  for(arg in names(formal.args)){ control[arg]<-list(get(arg)) }
  for(arg in names(p.args)){ control[arg]<-list(get(arg)) }

#suppressMessages(require(locfit, quietly=TRUE))
#require(coda)
out <- object$sample
if(is.null(out)){
  object$n <- min(object$x)
  object$lpriorm <- log(object$lprior)
}
if(!is.null(out)){
  outN <- out[,"N"]
  #a=locfit( ~ lp(outN, nn=0.35, h=0, maxk=500))
  a=locfit::locfit( ~ lp(outN,nn=0.5))
  xp <- seq(object$n,object$maxN, length=control$support)
  posdensN <- predict(a, newdata=xp)
  posdensN <- control$support*posdensN / ((object$maxN-object$n)*sum(posdensN))
  # Next from coda
  #hpd <- HPDinterval(object$sample[,"N"])[1:2]
  # MSH using locfit
  cy <- cumsum(posdensN/sum(posdensN))
  hpd <- c(xp[which.max(cy>((1-control$HPD.level)/2))],
           xp[which.max(cy>((1+control$HPD.level)/2))])
  if(is.na(hpd[1])) hpd[1] <- xp[1]
  if(is.na(hpd[2])) hpd[2] <- xp[length(xp)]
#
  map <- xp[which.max(posdensN)]
  mp <- sum(xp*posdensN)/sum(posdensN)
  l90 <- xp[which.max(cy>0.9)]
  l50 <- xp[which.max(cy>0.5)]
  l25 <- xp[which.max(cy>0.25)]
  l75 <- xp[which.max(cy>0.75)]
}
#
lpriorm <- exp(object$lpriorm-max(object$lpriorm))
lpriorm <- lpriorm[object$n+(1:length(lpriorm)) > object$n & object$n+(1:length(lpriorm)) < object$maxN]
lpriorm <- lpriorm / sum(lpriorm)
cy <- cumsum(lpriorm)
xp <- seq(object$n,object$maxN)
pl025 <- xp[which.max(cy>((1-control$HPD.level)/2))]
pl95  <- xp[which.max(cy>((1+control$HPD.level)/2))]
pl90  <- xp[which.max(cy>0.9)]
#
#cat(sprintf("Prior:\nMean = %d, Median = %d, Mode = %d, 90%% = %d, 25%% = %d, 75%% = %d.\n",
# round(object$mean.prior.size), round(object$median.prior.size), round(object$mode.prior.size), round(pl90), round(object$quartiles.prior.size[1]), round(object$quartiles.prior.size[2])))
#cat(sprintf("Posterior:\nMean = %d, Median = %d, MAP = %d, 90%% = %d, HPD = (%d, %d).\n",
# round(mp),round(l50),round(map),round(l90),round(hpd[1]),round(hpd[2])))
#
if(!is.null(out)){
 res <- matrix(c(
  round(object$mean.prior.size), round(object$median.prior.size), round(object$mode.prior.size),
  round(object$quartiles.prior.size[1]),round(object$quartiles.prior.size[2]),
  round(pl90),round(pl025),round(pl95),
  round(mp),round(l50),round(map),round(l25),round(l75),round(l90),round(hpd[1]),round(hpd[2])),byrow=TRUE,nrow=2)
  rownames(res) <- c("Prior","Posterior")
  colnames(res) <- c("Mean","Median","Mode","25%","75%","90%",
    paste(round(100*(1-control$HPD.level)/2,1),"%",sep=""),
    paste(round(100*(1+control$HPD.level)/2,1),"%",sep=""))
}else{
 res <- matrix(c(
  round(object$mean.prior.size), round(object$median.prior.size), round(object$mode.prior.size), round(object$quartiles.prior.size[1]), round(object$quartiles.prior.size[2]), round(pl90), round(pl025), round(pl95)
  ),byrow=TRUE,nrow=1)
  rownames(res) <- c("Prior")
  colnames(res) <- c("Mean","Median","Mode","25%","75%","90%",
    paste(round(100*(1-control$HPD.level)/2,1),"%",sep=""),
    paste(round(100*(1+control$HPD.level)/2,1),"%",sep=""))
}
res <- as.data.frame(res)
if(!is.null(out)){
  attr(res, "heading") <- "Summary of Population Size Estimation"
}else{
  attr(res, "heading") <- "Summary of Population Size Prior"
}
class(res) <- c("Anova","data.frame")
res
}
