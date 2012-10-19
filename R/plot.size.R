plot.size <- function(fit,data=NULL,xlim=2500){
require(locfit)
out <- fit$sample
outN <- out[,"N"]
a=locfit( ~ lp(outN, nn=2*0.35, h=0, maxk=500))
xp <- seq(fit$n,fit$maxN, by=1)
posdensN <- predict(a, newdata=xp)
posdensN <- posdensN / sum(posdensN)
lpriorm <- exp(fit$lpriorm-max(fit$lpriorm))
lpriorm <- lpriorm[fit$n+(1:length(lpriorm)) > fit$n & fit$n+(1:length(lpriorm)) < fit$maxN]
lpriorm <- lpriorm / sum(lpriorm)
#
plot(x=xp,y=posdensN,type='l', xlab="population size", 
  main="posterior for population size",
  ylim=c(0,max(posdensN,lpriorm)),
# sub="mean prior = 1000",
  ylab="posterior density",xlim=c(fit$n,xlim))
#
abline(v=fit$n,lty=2)
#
lpriorm <- exp(fit$lpriorm-max(fit$lpriorm))
lpriorm <- lpriorm/sum(lpriorm)
lines(x=fit$n+(1:length(lpriorm)),y=lpriorm,lty=2)
#
hpd <- HPDinterval(fit$sample[,"N"])[1:2]
#
abline(v=median(outN,na.rm=T),col=2)
abline(v=mean(outN,na.rm=T),col=3)
abline(v=c(fit$n,fit$maxN),lty=2)
abline(v=hpd,lty=2,col=4)
text(x=hpd[1],y=-0.000,col=4,cex=0.5,labels=paste(round(hpd[1])))
text(x=hpd[2],y=-0.000,col=4,cex=0.5,labels=paste(round(hpd[2])))
text(x=fit$n,y=0.000,labels=paste(fit$n),col=1,cex=0.5)
text(x=mean(outN,na.rm=T),y=-0.000,col=3,cex=0.5,labels=paste(round(mean(outN,na.rm=T))))
text(x=median(outN,na.rm=T),y=-0.000,col=2,cex=0.5,labels=paste(round(median(outN,na.rm=T))))
#
x <- fit$n+(1:length(lpriorm))
lines(x=fit$n+(1:length(lpriorm)),y=lpriorm,lty=2)
#
out[is.na(out)] <- apply(out,2,median,na.rm=T)
plot(density(out[,"mu"],na.rm=TRUE), xlab="mean network size", main="posterior for
mean network size in the population")
plot(density(out[,"sigma"],na.rm=TRUE), xlab="s.d. network size", main="posterior for s.d. of the network size")
#
if(!is.null(data)){
bbb <- tabulate(data,nbin=max(data))
bbb <- bbb/sum(bbb)
aaa <- barplot(bbb,names.arg=1:max(data),
  main="mean posterior network size distribution")
#abline(v=7)
lines(aaa, fit$predictive.data[1:max(data)], type='h', col='red', lwd=2)
}
invisible()
}
