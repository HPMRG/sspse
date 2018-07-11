#' Plot Summary and Diagnostics for Population Size Estimation Model Fits
#' 
#' This is the \code{plot} method for class \code{"sspse"}. Objects of
#' this class encapsulate the 
#' estimate of the posterior distribution of the
#' population size based on data collected by Respondent Driven Sampling. The
#' approach approximates the RDS via the Sequential Sampling model of Gile
#' (2008).  As such, it is referred to as the Sequential Sampling - Population Size Estimate (SS-PSE).
#' It uses the order of selection of the sample to provide information
#' on the distribution of network sizes over the population members.
#' 
#' By default it produces a density plot of the posterior for population size
#' and the prior for population size is overlaid. It also produces a
#' density plot of the posterior for mean network size in the population, the
#' posterior for standard deviation of the network size, and a density plot of
#' the posterior mean network size distribution with sample histogram overlaid.
#' 
#' @param x an object of class \code{"plot.sspse"}, usually, a result of a call
#' to \code{plot.sspse}.
#' @param xlim the (optional) x limits (x1, x2) of the plot of the posterior of
#' the population size.
#' @param data Optionally, the vector of degrees from the RDS in order they are
#' recorded and as passed to \code{\link{posteriorsize}}.
#' @param support the number of equally-spaced points to use for the support of
#' the estimated posterior density function.
#' @param HPD.level numeric; probability level of the highest probability
#' density interval determined from the estimated posterior.
#' @param N Optionally, an estimate of the population size to mark on the plots
#' as a reference point.
#' @param ylim the (optional) vertical limits (y1, y2) of the plot of the
#' posterior of the population size. A vertical axis is the probability density
#' scale.
#' @param mcmc logical; If TRUE, additionally create simple diagnostic plots
#' for the MCMC sampled statistics produced from the fit.
#' @param type character; This controls the types of plots produced.  If
#' \code{"N"}, a density plot of the posterior for population size is produced.
#' and the prior for population size is overlaid. If \code{"others"}, a 
#' density plot of the prior for population size, a
#' density plot of the posterior for mean network size in the population, the
#' posterior for standard deviation of the network size, and a density plot of
#' the posterior mean network size distribution with sample histogram overlaid
#' is produced.  If \code{"both"}, then all plots for \code{"N"} and
#' \code{"others"} are produced.
#' @param main an overall title for the posterior plot.
#' @param smooth the (optional) smoothing parameter for the density estimate.
#' @param \dots further arguments passed to or from other methods.
#' @seealso The model fitting function \code{\link{posteriorsize}},
#' \code{\link[graphics]{plot}}.
#' 
#' Function \code{\link[stats]{coef}} will extract the matrix of coefficients with
#' standard errors, t-statistics and p-values.
#' @keywords hplot
#' @references
#'
#' Gile, Krista J. (2008) \emph{Inference from Partially-Observed Network
#' Data}, Ph.D. Thesis, Department of Statistics, University of Washington.
#' 
#' Gile, Krista J. and Handcock, Mark S. (2010) \emph{Respondent-Driven
#' Sampling: An Assessment of Current Methodology}, Sociological Methodology
#' 40, 285-327.
#' 
#' Gile, Krista J. and Handcock, Mark S. (2014) \pkg{sspse}: Estimating Hidden 
#' Population Size using Respondent Driven Sampling Data
#' R package, Los Angeles, CA.  Version 0.5, \url{http://hpmrg.org}.
#' 
#' Handcock MS (2003).  \pkg{degreenet}: Models for Skewed Count Distributions
#' Relevant to Networks.  Statnet Project, Seattle, WA.  Version 1.2,
#' \url{http://statnetproject.org}.
#' 
#' Handcock, Mark S., Gile, Krista J. and Mar, Corinne M. (2014)
#' \emph{Estimating Hidden Population Size using Respondent-Driven Sampling
#' Data}, Electronic Journal of Statistics, 8, 1, 1491-1521
#' 
#' Handcock, Mark S., Gile, Krista J. and Mar, Corinne M. (2015)
#' \emph{Estimating the Size of Populations at High Risk for HIV using Respondent-Driven 
#' Sampling Data}, Biometrics.
#' @examples
#' 
#' N0 <- 200
#' n <- 100
#' K <- 10
#' 
#' # Create probabilities for a Waring distribution 
#' # with scaling parameter 3 and mean 5, but truncated at K=10.
#' probs <- c(0.33333333,0.19047619,0.11904762,0.07936508,0.05555556,
#'            0.04040404,0.03030303,0.02331002,0.01831502,0.01465201)
#' probs <- probs / sum(probs)
#' 
#' #
#' # Create a sample
#' #
#' set.seed(1)
#' pop<-sample(1:K, size=N0, replace = TRUE, prob = probs)
#' s<-sample(pop, size=n, replace = FALSE, prob = pop)
#'  
#' # Here interval=1 so that it will run faster. It should be higher in a 
#' # real application.
#' out <- posteriorsize(s=s,interval=1)
#' plot(out, HPD.level=0.9,data=pop[s])
#' summary(out, HPD.level=0.9)
#' # Let's look at some MCMC diagnostics
#' plot(out, HPD.level=0.9,mcmc=TRUE)
#' 
#' @importFrom coda mcmc
#' @method plot sspse
#' @export
plot.sspse <- function(x,
		       xlim=NULL,data=NULL,support=1000,HPD.level=0.90,N=NULL,ylim=NULL,mcmc=FALSE,type="both",
		       main="posterior for population size",smooth=4,...){
  p.args <- as.list( sys.call() )[-c(1,2)]
  formal.args<-formals(sys.function())[-c(1)]

  control<-list()
  names.formal.args <- names(formal.args)
  names.formal.args <- names.formal.args[-match("...",names.formal.args)]
  for(arg in names.formal.args){ control[arg]<-list(get(arg)) }
  for(arg in names(p.args)){ control[arg]<-list(get(arg)) }

out <- x$sample
# sabline <- function(v,...){graphics::segments(x0=v,...,y0=control$ylim[1],y1=control$ylim[2])}
sabline <- function(v,...){graphics::segments(x0=v,...,y0=graphics::par("usr")[3],y1=graphics::par("usr")[4])}
if(!is.null(out) & control$mcmc){
  mcmc.len <- min(1000, nrow(out))
  a=round(seq.int(from=1,to=nrow(out),length.out=mcmc.len))
  mcp <- attr(out,"mcpar")
  b=coda::mcmc(out[a,],start=mcp[1],end=mcp[2],thin=floor((mcp[2]-mcp[1])/mcmc.len + 1))
  graphics::plot(b)
  return(invisible())
}
#suppressMessages(require(locfit,quietly=TRUE))
#if(ask){graphics::par(ask=TRUE)}
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
  ##a=locfit( ~ lp(outN, nn=0.35, h=0, maxk=500))
  xp <- seq(x$n,x$maxN, length=control$support)
  posdensN=bgk_kde(data=outN,n=2^(ceiling(log(x$maxN-x$n)/log(2))),MIN=x$n,MAX=x$maxN, smooth=smooth)
  maxposdensN <- max(posdensN[1,],na.rm=TRUE)
  posdensN <- stats::spline(x=posdensN[1,],y=posdensN[2,],xout=xp)$y
  posdensN[xp > maxposdensN] <- 0
  #a=locfit::locfit( ~ lp(outN,nn=0.5))
  #posdensN <- predict(a, newdata=xp)
  posdensN <- control$support*posdensN / ((x$maxN-x$n)*sum(posdensN))
  #
  if(is.null(control$xlim)){control$xlim <- stats::quantile(outN,0.99)}
  if(is.null(control$ylim)){control$ylim <- c(0,max(posdensN,lpriorm))}
  if(control$type %in% c("N","both")){
  outp <- graphics::plot(x=xp,y=posdensN,type='l', xlab="population size", 
  # ylim=c(0,max(posdensN,lpriorm)),
  # sub="mean prior = 1000",
    ylab="posterior density",xlim=c(x$n,control$xlim),ylim=control$ylim, main=main)
  #
  graphics::legend('topright',lty=c(1,2,1,1,2,1),col=c(1,1,2,3,4,5),
    legend=c("posterior","prior","median","mean",
             paste(round(control$HPD.level*100),"% interval",sep=""),
             "mode"),
    bty="n",cex=0.75)
  sabline(v=x$n,lty=2)
  #
  lpriorm <- exp(x$lpriorm-max(x$lpriorm))
  lpriorm <- lpriorm/sum(lpriorm)
  graphics::lines(x=x$n+(1:length(lpriorm))-1,y=lpriorm,lty=2)
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
# sabline(v=stats::median(outN,na.rm=TRUE),col=2)
# sabline(v=mean(outN,na.rm=TRUE),col=3)
  sabline(v=mp,col=3)
  sabline(v=l50,col=2)
  sabline(v=map,col=5)
  sabline(v=c(x$n,x$maxN),lty=2)
  if(!is.null(control$N)){sabline(v=control$N,lty=1,col=1)}
  sabline(v=hpd,lty=2,col=4)
  yloc <- control$ylim[2]*0.06
  graphics::par(xpd=NA)
  graphics::text(x=hpd[1],y=-yloc,col=4,cex=0.5,labels=paste(round(hpd[1])))
  graphics::text(x=hpd[2],y=-yloc,col=4,cex=0.5,labels=paste(round(hpd[2])))
  graphics::text(x=x$n,y=-0,labels=paste(x$n),col=1,cex=1.0)
  graphics::text(x=round(mp),y=-yloc,col=3,cex=0.5,labels=paste(round(mp)))
  graphics::text(x=l50,y=-yloc,col=2,cex=0.5,labels=paste(round(l50)))
  graphics::text(x=map,y=-yloc,col=5,cex=0.5,labels=paste(round(map)))
  if(!is.null(control$N)){graphics::text(x=control$N,y=-yloc,col=1,cex=0.5,labels="truth")}
#
#cat(sprintf("Prior:\nMean = %d, Median = %d, Mode = %d, 25%% = %d, 75%% = %d.\n",
# round(x$mean.prior.size), round(x$median.prior.size), round(x$mode.prior.size), round(x$quartiles.prior.size[1]), round(x$quartiles.prior.size[2])))
#cat(sprintf("Posterior:\nMean = %d, Median = %d, MAP = %d, 90%% = %d, HPD = (%d, %d).\n",
# round(mp),round(l50),round(map),round(l90),round(hpd[1]),round(hpd[2])))
#
  graphics::lines(x=x$n+(1:length(lpriorm))-1,y=lpriorm,lty=2)
}
#
if(control$type %in% c("others","both")){
out[is.na(out)] <- apply(out,2,stats::median,na.rm=TRUE)
if("mu" %in% colnames(out)){
 graphics::plot(stats::density(out[,"mu"],na.rm=TRUE), xlab="mean network size", main="posterior for mean network size in the population")
 graphics::plot(stats::density(out[,"sigma"],na.rm=TRUE), xlab="s.d. network size", main="posterior for s.d. of the network size")
}
if("mu0" %in% colnames(out)){
 graphics::plot(stats::density(out[,"mu0"],na.rm=TRUE), xlab="mean network size for uninfected", main="posterior for mean network size in the uninfected population",sub="infected is dashed")
 graphics::lines(stats::density(out[,"mu1"],na.rm=TRUE),lty=2)
 graphics::plot(stats::density(out[,"sigma0"],na.rm=TRUE), xlab="s.d. network size", main="posterior for s.d. of the network size for uninfected",sub="infected is dashed")
 graphics::lines(stats::density(out[,"sigma1"],na.rm=TRUE),lty=2)
}
#
graphics::plot(seq_along(x$predictive.degree),y=x$predictive.degree, type='h',
col='red', lwd=2, xlab="degree",ylab="probability",
  main="mean posterior network size distribution")
if(!is.null(control$data)){
  Kmax <- max(seq_along(x$predictive.degree))
  bbb <- tabulate(control$data,nbins=Kmax) #, nbins=max(control$data))
  bbb <- bbb/sum(bbb)
  aaa <- graphics::barplot(bbb,names.arg=1:Kmax,add=FALSE,axes=TRUE,width=rep(0.5,length(bbb)),space=1,col=0,
    xlab="degree",ylab="probability", xlim=c(1,Kmax),
    main="posterior with sample histogram overlaid")
  graphics::lines(x=-0.25+seq_along(x$predictive.degree),y=x$predictive.degree, type='h', col='red', lwd=2)
}}
}else{
if(control$type %in% c("others","both")){
  cy <- cumsum(lpriorm)
  #
  if(is.null(control$xlim)){control$xlim <- xp[which.max(cy>0.99)]}
  if(is.null(control$ylim)){control$ylim <- c(0,max(lpriorm))}
  graphics::plot(x=xp,y=lpriorm,type='l', xlab="population size", 
    main="prior for population size",
    ylab="prior density",xlim=c(x$n,control$xlim),ylim=control$ylim)
  #
  sabline(v=x$n,lty=2)
  #
#
  map <- xp[which.max(lpriorm)]
  mp <- sum(xp*lpriorm)/sum(lpriorm)
  l90 <- xp[which.max(cy>0.9)]
  l50 <- xp[which.max(cy>0.5)]
  hpd <- c(xp[which.max(cy>0.25)],xp[which.max(cy>0.75)])
#
  sabline(v=l50,col=2)
  sabline(v=mp,col=3)
  sabline(v=map,col=5)
  sabline(v=c(x$n,x$maxN),lty=2)
  if(!is.null(control$N)){sabline(v=control$N,lty=1,col=1)}
  sabline(v=hpd,lty=2,col=4)
  yloc <- control$ylim[2]*0.06
  graphics::par(xpd=NA)
  graphics::text(x=hpd[1],y=-yloc,col=4,cex=0.5,labels=paste(round(hpd[1])))
  graphics::text(x=hpd[2],y=-yloc,col=4,cex=0.5,labels=paste(round(hpd[2])))
  graphics::text(x=x$n,y=-0.0,labels=paste(x$n),col=1,cex=1.0)
  graphics::text(x=mp,y=-yloc,col=3,cex=0.5,labels=paste(round(mp)))
  graphics::text(x=l50,y=-yloc,col=2,cex=0.5,labels=paste(round(l50)))
  graphics::text(x=map,y=-yloc,col=5,cex=0.5,labels=paste(round(map)))
  if(!is.null(control$N)){graphics::text(x=control$N,y=-yloc,col=1,cex=0.5,labels="truth")}
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
