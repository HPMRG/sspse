summary.size <- function(fit,support=1000,HPD.level=0.95){
require(locfit)
#require(coda)
out <- fit$sample
if(is.null(out)){
  fit$n <- min(fit$x)
  fit$lpriorm <- log(fit$lprior)
}
if(!is.null(out)){
  outN <- out[,"N"]
  #a=locfit( ~ lp(outN, nn=0.35, h=0, maxk=500))
  a=locfit( ~ lp(outN,nn=0.5))
  xp <- seq(fit$n,fit$maxN, length=support)
  posdensN <- predict(a, newdata=xp)
  posdensN <- support*posdensN / ((fit$maxN-fit$n)*sum(posdensN))
  # Next from coda
  #hpd <- HPDinterval(fit$sample[,"N"])[1:2]
  # MSH using locfit
  cy <- cumsum(posdensN/sum(posdensN))
  hpd <- c(xp[which.max(cy>((1-HPD.level)/2))],
           xp[which.max(cy>((1+HPD.level)/2))])
  if(is.na(hpd[1])) hpd[1] <- xp[1]
  if(is.na(hpd[2])) hpd[2] <- xp[length(xp)]
#
  map <- xp[which.max(posdensN)]
  mp <- sum(xp*posdensN)/sum(posdensN)
  l90 <- xp[which.max(cy>0.9)]
  l50 <- xp[which.max(cy>0.5)]
}
#
lpriorm <- exp(fit$lpriorm-max(fit$lpriorm))
lpriorm <- lpriorm[fit$n+(1:length(lpriorm)) > fit$n & fit$n+(1:length(lpriorm)) < fit$maxN]
lpriorm <- lpriorm / sum(lpriorm)
cy <- cumsum(lpriorm)
xp <- seq(fit$n,fit$maxN)
pl90 <- xp[which.max(cy>0.9)]
#
#cat(sprintf("Prior:\nMean = %d, Median = %d, Mode = %d, 90%% = %d, 25%% = %d, 75%% = %d.\n",
# round(fit$mean.prior.size), round(fit$median.prior.size), round(fit$mode.prior.size), round(pl90), round(fit$quartiles.prior.size[1]), round(fit$quartiles.prior.size[2])))
#cat(sprintf("Posterior:\nMean = %d, Median = %d, MAP = %d, 90%% = %d, HPD = (%d, %d).\n",
# round(mp),round(l50),round(map),round(l90),round(hpd[1]),round(hpd[2])))
#
if(!is.null(out)){
 res <- matrix(c(
  round(fit$mean.prior.size), round(fit$median.prior.size), round(fit$mode.prior.size), round(pl90), round(fit$quartiles.prior.size[1]), round(fit$quartiles.prior.size[2]),
  round(mp),round(l50),round(map),round(l90),round(hpd[1]),round(hpd[2])),byrow=TRUE,nrow=2)
  rownames(res) <- c("Prior","Posterior")
}else{
 res <- matrix(c(
  round(fit$mean.prior.size), round(fit$median.prior.size), round(fit$mode.prior.size), round(pl90), round(fit$quartiles.prior.size[1]), round(fit$quartiles.prior.size[2])
  ),byrow=TRUE,nrow=1)
  rownames(res) <- c("Prior")
}
colnames(res) <- c("Mean","Median","Mode","90%","25%","75%")
res <- as.data.frame(res)
if(!is.null(out)){
  attr(res, "heading") <- "Summary of Population Size Estimation"
}else{
  attr(res, "heading") <- "Summary of Population Size Prior"
}
print.Anova(res)
invisible(res)
}
