plot.psess <- function(x,
		       xlim=NULL,data=NULL,support=1000,HPD.level=0.90,N=NULL,ylim=NULL,mcmc=FALSE,type="both",...){
  p.args <- as.list( sys.call() )[-c(1,2)]
  formal.args<-formals(sys.function())[-c(1)]

  control<-list()
  names.formal.args <- names(formal.args)
  names.formal.args <- names.formal.args[-match("...",names.formal.args)]
  for(arg in names.formal.args){ control[arg]<-list(get(arg)) }
  for(arg in names(p.args)){ control[arg]<-list(get(arg)) }

out <- x$sample
if(!is.null(out) & control$mcmc){
  mcmc.len <- min(1000, nrow(out))
  a=round(seq.int(from=1,to=nrow(out),length.out=mcmc.len))
  mcp <- attr(out,"mcpar")
  b=coda::mcmc(out[a,],start=mcp[1],end=mcp[2],thin=floor((mcp[2]-mcp[1])/mcmc.len + 1))
  plot(b)
  return(invisible())
}
#suppressMessages(require(locfit,quietly=TRUE))
#if(ask){par(ask=TRUE)}
if(is.null(out)){
  x$n <- min(x$x)
  x$lpriorm <- log(x$lprior)
}
xp <- x$n+(1:length(x$lpriorm))-1
lpriorm <- exp(x$lpriorm-max(x$lpriorm))
lpriorm <- lpriorm[xp >= x$n & xp <= x$maxN]
lpriorm <- lpriorm / sum(lpriorm)
xp <- x$n+(1:length(lpriorm))-1
if(!is.null(out)){
  outN <- out[,"N"]
  #a=locfit( ~ lp(outN, nn=0.35, h=0, maxk=500))
  a=locfit::locfit( ~ lp(outN,nn=0.5))
  xp <- seq(x$n,x$maxN, length=control$support)
  posdensN <- predict(a, newdata=xp)
  posdensN <- control$support*posdensN / ((x$maxN-x$n)*sum(posdensN))
  #
  if(is.null(control$xlim)){control$xlim <- quantile(outN,0.99)}
  if(is.null(control$ylim)){control$ylim <- c(0,max(posdensN,lpriorm))}
  if(control$type %in% c("N","both")){
  plot(x=xp,y=posdensN,type='l', xlab="population size", 
    main="posterior for population size",
  # ylim=c(0,max(posdensN,lpriorm)),
  # sub="mean prior = 1000",
    ylab="posterior density",xlim=c(x$n,control$xlim),ylim=control$ylim)
  #
  legend('topright',lty=c(1,2,1,1,1,1,1),col=c(1,1,2,3,4,5),
    legend=c("posterior","prior","median","mean",
             paste(round(control$HPD.level*100),"% interval",sep=""),
             "mode"),
    bty="n",cex=0.75)
  abline(v=x$n,lty=2)
  #
  lpriorm <- exp(x$lpriorm-max(x$lpriorm))
  lpriorm <- lpriorm/sum(lpriorm)
  lines(x=x$n+(1:length(lpriorm))-1,y=lpriorm,lty=2)
  # Next from coda
  #hpd <- HPDinterval(x$sample[,"N"])[1:2]
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
#
# abline(v=median(outN,na.rm=TRUE),col=2)
# abline(v=mean(outN,na.rm=TRUE),col=3)
  abline(v=l50,col=2)
  abline(v=mp,col=3)
  abline(v=c(x$n,x$maxN),lty=2)
  if(!is.null(control$N)){abline(v=control$N,lty=1,col=1)}
  abline(v=hpd,lty=2,col=4)
  yloc <- control$ylim[2]*0.06
  par(xpd=NA)
  text(x=hpd[1],y=-yloc,col=4,cex=0.5,labels=paste(round(hpd[1])))
  text(x=hpd[2],y=-yloc,col=4,cex=0.5,labels=paste(round(hpd[2])))
  text(x=x$n,y=-0,labels=paste(x$n),col=1,cex=1.0)
  text(x=round(mp),y=-yloc,col=3,cex=0.5,labels=paste(round(mp)))
  text(x=l50,y=-yloc,col=2,cex=0.5,labels=paste(round(l50)))
  text(x=map,y=-yloc,col=5,cex=0.5,labels=paste(round(map)))
  if(!is.null(control$N)){text(x=control$N,y=-yloc,col=1,cex=0.5,labels="truth")}
#
#cat(sprintf("Prior:\nMean = %d, Median = %d, Mode = %d, 25%% = %d, 75%% = %d.\n",
# round(x$mean.prior.size), round(x$median.prior.size), round(x$mode.prior.size), round(x$quartiles.prior.size[1]), round(x$quartiles.prior.size[2])))
#cat(sprintf("Posterior:\nMean = %d, Median = %d, MAP = %d, 90%% = %d, HPD = (%d, %d).\n",
# round(mp),round(l50),round(map),round(l90),round(hpd[1]),round(hpd[2])))
#
lines(x=x$n+(1:length(lpriorm))-1,y=lpriorm,lty=2)
}
#
if(control$type %in% c("others","both")){
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
plot(seq_along(x$predictive.degree),y=x$predictive.degree, type='h',
col='red', lwd=2, xlab="degree",ylab="probability",
  main="mean posterior network size distribution")
if(!is.null(control$data)){
  Kmax <- max(seq_along(x$predictive.degree))
  bbb <- tabulate(control$data,nbins=Kmax) #, nbins=max(control$data))
  bbb <- bbb/sum(bbb)
  aaa <- barplot(bbb,names.arg=1:Kmax,add=FALSE,axes=TRUE,width=rep(0.5,length(bbb)),space=1,col=0,
    xlab="degree",ylab="probability", xlim=c(1,Kmax),
    main="posterior with sample histogram overlaid")
  lines(x=-0.25+seq_along(x$predictive.degree),y=x$predictive.degree, type='h', col='red', lwd=2)
}}
}else{
if(control$type %in% c("others","both")){
  cy <- cumsum(lpriorm)
  #
  if(is.null(control$xlim)){control$xlim <- xp[which.max(cy>0.99)]}
  if(is.null(control$ylim)){control$ylim <- c(0,max(lpriorm))}
  plot(x=xp,y=lpriorm,type='l', xlab="population size", 
    main="prior for population size",
    ylab="prior density",xlim=c(x$n,control$xlim),ylim=control$ylim)
  #
  abline(v=x$n,lty=2)
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
  abline(v=c(x$n,x$maxN),lty=2)
  if(!is.null(control$N)){abline(v=control$N,lty=1,col=1)}
  abline(v=hpd,lty=2,col=4)
  yloc <- control$ylim[2]*0.06
  par(xpd=NA)
  text(x=hpd[1],y=-yloc,col=4,cex=0.5,labels=paste(round(hpd[1])))
  text(x=hpd[2],y=-yloc,col=4,cex=0.5,labels=paste(round(hpd[2])))
  text(x=x$n,y=-0.0,labels=paste(x$n),col=1,cex=1.0)
  text(x=mp,y=-yloc,col=3,cex=0.5,labels=paste(round(mp)))
  text(x=l50,y=-yloc,col=2,cex=0.5,labels=paste(round(l50)))
  text(x=map,y=-yloc,col=5,cex=0.5,labels=paste(round(map)))
  if(!is.null(control$N)){text(x=control$N,y=-yloc,col=1,cex=0.5,labels="truth")}
}
#
#cat(sprintf("Prior:\nMean = %d, Median = %d, Mode = %d, 25%% = %d, 75%% = %d.\n",
# round(x$mean.prior.size), round(x$median.prior.size), round(x$mode.prior.size), round(x$quartiles.prior.size[1]), round(x$quartiles.prior.size[2])))
#cat(sprintf("Posterior:\nMean = %d, Median = %d, MAP = %d, 90%% = %d, HPD = (%d, %d).\n",
# round(mp),round(l50),round(map),round(l90),round(hpd[1]),round(hpd[2])))
#
}
invisible()
}
