plot.size <-
function(fit,xlim=NULL,data=NULL,support=1000,HPD.level=0.95,N=NULL,ylim=NULL){
require(locfit)
require(coda)
out <- fit$sample
outN <- out[,"N"]
#a=locfit( ~ lp(outN, nn=0.35, h=0, maxk=500))
a=locfit( ~ lp(outN,nn=0.5))
xp <- seq(fit$n,fit$maxN, length=support)
posdensN <- predict(a, newdata=xp)
posdensN <- support*posdensN / ((fit$maxN-fit$n)*sum(posdensN))
lpriorm <- exp(fit$lpriorm-max(fit$lpriorm))
lpriorm <- lpriorm[fit$n+(1:length(lpriorm)) > fit$n & fit$n+(1:length(lpriorm)) < fit$maxN]
lpriorm <- lpriorm / sum(lpriorm)
#
if(is.null(xlim)){xlim <- quantile(outN,0.99)}
if(is.null(ylim)){ylim <- c(0,max(posdensN,lpriorm))}
plot(x=xp,y=posdensN,type='l', xlab="population size", 
  main="posterior for population size",
# ylim=c(0,max(posdensN,lpriorm)),
# sub="mean prior = 1000",
  ylab="posterior density",xlim=c(fit$n,xlim),ylim=c(0,ylim))
#
abline(v=fit$n,lty=2)
#
lpriorm <- exp(fit$lpriorm-max(fit$lpriorm))
lpriorm <- lpriorm/sum(lpriorm)
lines(x=fit$n+(1:length(lpriorm)),y=lpriorm,lty=2)
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
#
abline(v=median(outN,na.rm=T),col=2)
abline(v=mean(outN,na.rm=T),col=3)
abline(v=c(fit$n,fit$maxN),lty=2)
if(!is.null(N)){abline(v=N,lty=1,col=1)}
abline(v=hpd,lty=2,col=4)
text(x=hpd[1],y=-0.000,col=4,cex=0.5,labels=paste(round(hpd[1])))
text(x=hpd[2],y=-0.000,col=4,cex=0.5,labels=paste(round(hpd[2])))
text(x=fit$n,y=0.000,labels=paste(fit$n),col=1,cex=0.5)
text(x=mean(outN,na.rm=T),y=-0.000,col=3,cex=0.5,labels=paste(round(mean(outN,na.rm=T))))
text(x=median(outN,na.rm=T),y=-0.000,col=2,cex=0.5,labels=paste(round(median(outN,na.rm=T))))
text(x=map,y=-0.000,col=5,cex=0.5,labels=paste(round(map)))
if(!is.null(N)){text(x=N,y=-0.000,col=1,cex=0.5,labels="truth")}
#
cat(sprintf("Prior:\nMean = %d, Median = %d, Mode = %d, 25%% = %d, 75%% = %d.\n",
 round(fit$mean.prior.size), round(fit$median.prior.size), round(fit$mode.prior.size), round(fit$quartiles.prior.size[1]), round(fit$quartiles.prior.size[2])))
cat(sprintf("Posterior:\nMean = %d, Median = %d, MAP = %d, 90%% = %d, HPD = (%d, %d).\n",
 round(mp),round(l50),round(map),round(l90),round(hpd[1]),round(hpd[2])))
#
x <- fit$n+(1:length(lpriorm))
lines(x=fit$n+(1:length(lpriorm)),y=lpriorm,lty=2)
#
out[is.na(out)] <- apply(out,2,median,na.rm=T)
plot(density(out[,"mu"],na.rm=TRUE), xlab="mean network size", main="posterior for mean network size in the population")
plot(density(out[,"sigma"],na.rm=TRUE), xlab="s.d. network size", main="posterior for s.d. of the network size")
#
plot(seq_along(fit$predictive.degree),y=fit$predictive.degree, type='h',
col='red', lwd=2, xlab="degree",ylab="probability",
  main="mean posterior network size distribution")
if(!is.null(data)){
  bbb <- tabulate(data,nbins=max(data))
  bbb <- bbb/sum(bbb)
  aaa <- barplot(bbb,names.arg=1:max(data),add=T)
}
invisible()
}
