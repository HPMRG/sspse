plot.size <-
function(fit,xlim=NULL,data=NULL,support=1000,HPD.level=0.95,N=NULL,ylim=NULL,mcmc=FALSE,type="both"){
out <- fit$sample
if(!is.null(out) & mcmc){
# suppressPackageStartupMessages(require(coda, quietly=TRUE))
  mcmc.len <- min(1000, nrow(out))
  a=round(seq.int(from=1,to=nrow(out),length=mcmc.len))
  mcp <- attr(out,"mcpar")
  b=coda::mcmc(out[a,],start=mcp[1],end=mcp[2],thin=round((mcp[2]-mcp[1])/mcmc.len))
  plot(b)
  return(invisible())
}
suppressMessages(require(locfit,quietly=TRUE))
#if(ask){par(ask=TRUE)}
if(is.null(out)){
  fit$n <- min(fit$x)
  fit$lpriorm <- log(fit$lprior)
}
xp <- fit$n+(1:length(fit$lpriorm))-1
lpriorm <- exp(fit$lpriorm-max(fit$lpriorm))
lpriorm <- lpriorm[xp >= fit$n & xp <= fit$maxN]
lpriorm <- lpriorm / sum(lpriorm)
xp <- fit$n+(1:length(lpriorm))-1
x <- fit$n+(1:length(lpriorm))-1
if(!is.null(out)){
  outN <- out[,"N"]
  #a=locfit( ~ lp(outN, nn=0.35, h=0, maxk=500))
  a=locfit::locfit( ~ lp(outN,nn=0.5))
  xp <- seq(fit$n,fit$maxN, length=support)
  posdensN <- predict(a, newdata=xp)
  posdensN <- support*posdensN / ((fit$maxN-fit$n)*sum(posdensN))
  #
  if(is.null(xlim)){xlim <- quantile(outN,0.99)}
  if(is.null(ylim)){ylim <- c(0,max(posdensN,lpriorm))}
  if(type %in% c("N","both")){
  plot(x=xp,y=posdensN,type='l', xlab="population size", 
    main="posterior for population size",
  # ylim=c(0,max(posdensN,lpriorm)),
  # sub="mean prior = 1000",
    ylab="posterior density",xlim=c(fit$n,xlim),ylim=ylim)
  #
  abline(v=fit$n,lty=2)
  #
  lpriorm <- exp(fit$lpriorm-max(fit$lpriorm))
  lpriorm <- lpriorm/sum(lpriorm)
  lines(x=fit$n+(1:length(lpriorm))-1,y=lpriorm,lty=2)
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
  abline(v=median(outN,na.rm=TRUE),col=2)
  abline(v=mean(outN,na.rm=TRUE),col=3)
  abline(v=c(fit$n,fit$maxN),lty=2)
  if(!is.null(N)){abline(v=N,lty=1,col=1)}
  abline(v=hpd,lty=2,col=4)
  text(x=hpd[1],y=-0.000,col=4,cex=1.0,labels=paste(round(hpd[1])))
  text(x=hpd[2],y=-0.000,col=4,cex=1.0,labels=paste(round(hpd[2])))
  text(x=fit$n,y=0.000,labels=paste(fit$n),col=1,cex=1.0)
  text(x=mean(outN,na.rm=TRUE),y=-0.000,col=3,cex=1.0,labels=paste(round(mean(outN,na.rm=TRUE))))
  text(x=median(outN,na.rm=TRUE),y=-0.000,col=2,cex=1.0,labels=paste(round(median(outN,na.rm=TRUE))))
  text(x=map,y=-0.000,col=5,cex=1.0,labels=paste(round(map)))
  if(!is.null(N)){text(x=N,y=-0.000,col=1,cex=1.0,labels="truth")}
#
#cat(sprintf("Prior:\nMean = %d, Median = %d, Mode = %d, 25%% = %d, 75%% = %d.\n",
# round(fit$mean.prior.size), round(fit$median.prior.size), round(fit$mode.prior.size), round(fit$quartiles.prior.size[1]), round(fit$quartiles.prior.size[2])))
#cat(sprintf("Posterior:\nMean = %d, Median = %d, MAP = %d, 90%% = %d, HPD = (%d, %d).\n",
# round(mp),round(l50),round(map),round(l90),round(hpd[1]),round(hpd[2])))
#
lines(x=fit$n+(1:length(lpriorm))-1,y=lpriorm,lty=2)
}
#
if(type %in% c("others","both")){
out[is.na(out)] <- apply(out,2,median,na.rm=TRUE)
if("mu" %in% colnames(out)){
 plot(density(out[,"mu"],na.rm=TRUE), xlab="mean network size", main="posterior for mean network size in the population")
 plot(density(out[,"sigma"],na.rm=TRUE), xlab="s.d. network size", main="posterior for s.d. of the network size")
}
if("mu0" %in% colnames(out)){
 plot(density(out[,"mu0"],na.rm=TRUE), xlab="mean network size for uninfected", main="posterior for mean network size in the uninfected population",sub="infected is dashed")
 lines(density(out[,"mu1"],na.rm=TRUE),lty=2)
 plot(density(out[,"sigma0"],na.rm=TRUE), xlab="s.d. network size", main="posterior for s.d. of the network size for uninfected",sub="infected is dashed")
 lines(density(out[,"sigma1"],na.rm=TRUE),lty=2)
}
#
plot(seq_along(fit$predictive.degree),y=fit$predictive.degree, type='h',
col='red', lwd=2, xlab="degree",ylab="probability",
  main="mean posterior network size distribution")
if(!is.null(data)){
  Kmax <- max(seq_along(fit$predictive.degree))
  bbb <- tabulate(data,nbins=Kmax) #, nbins=max(data))
  bbb <- bbb/sum(bbb)
  aaa <- barplot(bbb,names.arg=1:Kmax,add=FALSE,axes=TRUE,width=rep(0.5,length(bbb)),space=1,col=0,
    xlab="degree",ylab="probability", xlim=c(1,Kmax),
    main="posterior with sample histogram overlaid")
  lines(x=-0.25+seq_along(fit$predictive.degree),y=fit$predictive.degree, type='h', col='red', lwd=2)
}}
}else{
if(type %in% c("others","both")){
  cy <- cumsum(lpriorm)
  #
  if(is.null(xlim)){xlim <- xp[which.max(cy>0.99)]}
  if(is.null(ylim)){ylim <- c(0,max(lpriorm))}
  plot(x=xp,y=lpriorm,type='l', xlab="population size", 
    main="prior for population size",
    ylab="prior density",xlim=c(fit$n,xlim),ylim=ylim)
  #
  abline(v=fit$n,lty=2)
  #
#
  map <- xp[which.max(lpriorm)]
  mp <- sum(xp*lpriorm)/sum(lpriorm)
  l90 <- xp[which.max(cy>0.9)]
  l50 <- xp[which.max(cy>0.5)]
  hpd <- c(xp[which.max(cy>0.25)],xp[which.max(cy>0.75)])
#
  abline(v=l50,col=2)
  abline(v=mp,col=3)
  abline(v=c(fit$n,fit$maxN),lty=2)
  if(!is.null(N)){abline(v=N,lty=1,col=1)}
  abline(v=hpd,lty=2,col=4)
  text(x=hpd[1],y=-0.000,col=4,cex=1.0,labels=paste(round(hpd[1])))
  text(x=hpd[2],y=-0.000,col=4,cex=1.0,labels=paste(round(hpd[2])))
  text(x=fit$n,y=0.000,labels=paste(fit$n),col=1,cex=1.0)
  text(x=mp,y=-0.000,col=3,cex=1.0,labels=paste(round(mp)))
  text(x=l50,y=-0.000,col=2,cex=1.0,labels=paste(round(l50)))
  text(x=map,y=-0.000,col=5,cex=1.0,labels=paste(round(map)))
  if(!is.null(N)){text(x=N,y=-0.000,col=1,cex=1.0,labels="truth")}
}
#
#cat(sprintf("Prior:\nMean = %d, Median = %d, Mode = %d, 25%% = %d, 75%% = %d.\n",
# round(fit$mean.prior.size), round(fit$median.prior.size), round(fit$mode.prior.size), round(fit$quartiles.prior.size[1]), round(fit$quartiles.prior.size[2])))
#cat(sprintf("Posterior:\nMean = %d, Median = %d, MAP = %d, 90%% = %d, HPD = (%d, %d).\n",
# round(mp),round(l50),round(map),round(l90),round(hpd[1]),round(hpd[2])))
#
}
invisible()
}
