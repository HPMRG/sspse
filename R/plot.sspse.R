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
#' and the prior for population size is overlaid. If \code{"summary"}, a 
#' density plot of the posterior for mean visibility in the population and
#' a plot of the posterior for standard deviation of the visibility in the population.
#' If \code{"visibility"}, a density plot of
#' the visibility distribution (its posterior mean) and the same plot with the 
#' with visibilities of those in the sample overlaid.
#' If \code{"degree"}, a scatter plot of the visibilities verses the reported network sizes for 
#' those in the sample.
#' If \code{"prior"}, a density plot of the prior for population size is produced.
#' If \code{"all"}, then all plots for \code{"N"}, \code{"summary"}, \code{"visibility"} and
#' \code{"degree"} are produced.
#' In all cases the visibilities are estimated (by their posterior means).
#' @param main an overall title for the posterior plot.
#' @param smooth the (optional) smoothing parameter for the density estimate.
#' @param include.tree logical; If \code{TRUE}, 
#' augment the reported network size by the number of recruits and one for the recruiter (if any).
#' This reflects a more accurate value for the visibility, but is not the reported degree.
#' In particular, it typically produces a positive visibility (compared to a possibility zero reported degree). 
#' @param cex.main an overall title for the posterior plot.
#' @param log.degree a character string which contains \code{"x"} if the (horizontal) degree axis in the plot
#' of the estimated visibilites for each respondent verses their reported network sizes be logarithmic. 
#' A value of \code{"y"} uses a logarithmic visibility axis and \code{"xy"} both. The default is \code{""}, no logarithmic axes.
#' @param method character; The method to use for density estimation (default Gaussian Kernel; "bgk").
#' "Bayes" uses a Bayesian density estimator which has good properties.
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
#' R package, Los Angeles, CA.  Version 0.5, \url{https://hpmrg.org}.
#' 
#' Handcock MS (2003).  \pkg{degreenet}: Models for Skewed Count Distributions
#' Relevant to Networks.  Statnet Project, Seattle, WA.  Version 1.2,
#' \url{https://statnet.org}.
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
#' \dontrun{
#' data(fauxmadrona)
#' # Here interval=1 and samplesize=50 so that it will run faster. It should be much higher
#' # in a real application.
#' fit <- posteriorsize(fauxmadrona, median.prior.size=1000,
#'                                   warmup=10, interval=1, samplesize=50)
#' summary(fit)
#' # Let's look at some MCMC diagnostics
#' plot(fit, mcmc=TRUE)
#' }
#' 
#' @importFrom coda mcmc
#' @method plot sspse
#' @export
plot.sspse <- function(x,
		       xlim=NULL,support=1000,HPD.level=0.90,N=NULL,ylim=NULL,mcmc=FALSE,type="all",
		       main="Posterior for population size",smooth=4,include.tree=TRUE,cex.main=1,log.degree="",method="bgk",...){
  p.args <- as.list( sys.call() )[-c(1,2)]
  formal.args<-formals(sys.function())[-c(1)]

 #control<-list()
  control<-list(samples=4000, warmup=1000)
  names.formal.args <- names(formal.args)
  names.formal.args <- names.formal.args[-match("...",names.formal.args)]
  for(arg in names.formal.args){ control[arg]<-list(get(arg)) }
  for(arg in names(p.args)){ control[arg]<-list(get(arg)) }

  if(control$mcmc){control$type <- "mcmc"}

out <- x$sample
# Remove NaN and NA by replacing with the minimum value
out <- apply(out,2,function(x){x[is.na(x)] <- min(x,na.rm=TRUE);x})
attr(out,"mcpar") <- attr(x$sample,"mcpar")
attr(out,"class") <- attr(x$sample,"class")
# sabline <- function(v,...){graphics::segments(x0=v,...,y0=control$ylim[1],y1=control$ylim[2])}
sabline <- function(v,...){graphics::segments(x0=v,...,y0=graphics::par("usr")[3],y1=graphics::par("usr")[4])}
if(!is.null(out) & control$type == "mcmc"){
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
  if(method=="bgk" || !requireNamespace("densEstBayes", quietly = TRUE)){
#   posdensN=bgk_kde(data=outN,n=2^(ceiling(log(x$maxN-x$n)/log(2))),MIN=x$n,MAX=x$maxN, smooth=smooth)
#   maxposdensN <- max(posdensN[1,],na.rm=TRUE)
#   posdensN <- stats::spline(x=posdensN[1,],y=posdensN[2,],xout=xp)$y
    posdensN=KernSmooth::bkde(x=log(outN), kernel = "normal", gridsize = length(xp), range.x=log(c(x$n,x$maxN)))
#   posdensN=density(x=log(outN),n=control$support,from=log(x$n),to=log(x$maxN))
#   xp=exp(posdensN$x)
#   posdensN=posdensN$y/xp
    maxposdensN <- max(exp(posdensN$x),na.rm=TRUE)
    posdensN <- stats::spline(x=exp(posdensN$x),y=posdensN$y/exp(posdensN$x),xout=xp)$y
    posdensN[xp > maxposdensN] <- 0
    #a=locfit::locfit( ~ lp(outN,nn=0.5))
    #posdensN <- predict(a, newdata=xp)
  }else{
    a=densEstBayes::densEstBayes(outN,method="NUTS",
      control=densEstBayes::densEstBayes.control(range.x=c(x$n*0.95,x$maxN*1.05),#numBins=min(401,round(length(outN)/2)),
                                                 nKept=control$samples,nWarm=control$warmup))
    xTrang <- seq(-0.05, 1.05, length = length(xp))
    Xg <- cbind(1,xTrang)
    Zg <- .ZOSull(xTrang,intKnots=a$intKnots,range.x=c(-0.05,1.05))
    Cg <- cbind(Xg,Zg)
    betauMCMC <- a$stochaFitObj$betauMCMC
    etaHatMCMC <- crossprod(t(Cg),betauMCMC)
    posdensN <- exp(apply(etaHatMCMC, 1, mean))
  }
  posdensN <- control$support*posdensN / ((x$maxN-x$n)*sum(posdensN))
  #
  if(is.null(control$xlim)){control$xlim <- stats::quantile(outN,0.99)}
  if(is.null(control$ylim)){control$ylim <- c(0,max(posdensN,lpriorm))}
  if(control$type %in% c("N","all")){
  outp <- graphics::plot(x=xp,y=posdensN,type='l', xlab="population size", 
  # ylim=c(0,max(posdensN,lpriorm)),
  # sub="mean prior = 1000",
    ylab="Posterior Density",xlim=c(x$n,control$xlim),ylim=control$ylim, main=main, cex.main=cex.main)
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
if(control$type %in% c("summary","all")){
 out[is.na(out)] <- apply(out,2,stats::median,na.rm=TRUE)
 if("mu" %in% colnames(out)){
  graphics::plot(stats::density(out[,"mu"],na.rm=TRUE), xlab="mean visibility", main="Posterior for mean visibility in the population", cex.main=cex.main)
  graphics::plot(stats::density(out[,"sigma"],na.rm=TRUE), xlab="s.d. visibility", main="Posterior for s.d. of the visibility", cex.main=cex.main)
 }
 if("mu0" %in% colnames(out)){ # posteriordisease
  graphics::plot(stats::density(out[,"mu0"],na.rm=TRUE), xlab="mean visibility for uninfected", xlim=range(c(out[,"mu0"],out[,"mu1"])),
    main="Posterior for mean visibility in the uninfected population",sub="infected is dashed", cex.main=cex.main)
  graphics::lines(stats::density(out[,"mu1"],na.rm=TRUE),col=3)
  graphics::legend('topleft',col=c(1,3), lty=c(1,1),
    legend=c("uninfected","infected"), bty="n",cex=0.75)
  graphics::plot(stats::density(out[,"sigma0"],na.rm=TRUE), xlab="s.d. visibility", xlim=range(c(out[,"sigma0"],out[,"sigma1"])),
    main="Posterior for s.d. of the visibility for uninfected",sub="infected is dashed", cex.main=cex.main)
  graphics::lines(stats::density(out[,"sigma1"],na.rm=TRUE),col=3)
  graphics::legend('topleft',col=c(1,3), lty=c(1,1),
    legend=c("uninfected","infected"), bty="n",cex=0.75)
 }
}
#
if(control$type %in% c("visibility","degree","all")){
 if(control$type %in% c("visibility")){
  graphics::plot(seq_along(x$predictive.visibility),y=x$predictive.visibility, type='h',
  col='red', lwd=2, xlab="visibility",ylab="probability", ylim=c(0,max(x$predictive.visibility)),
    main="Visibility distribution", cex.main=cex.main)
 }
 if(!is.null(x$data)){
   if(methods::is(x$data,"rds.data.frame")){
    if(is.null(attr(x$data,"network.size.variable"))){
      stop("Passed data must have a network.size attribute.")
    }
#   nw <- RDS::get.wave(x$data)
#   ns <- RDS::get.seed.id(rds.datax$data

    network.size <- as.numeric(x$data[[attr(x$data,"network.size.variable")]])
    #Augment the reported network size by the number of recruits and the recruiter (if any).
    if(include.tree){
     nr <- RDS::get.number.of.recruits(x$data)
     is.seed <- (RDS::get.rid(x$data)=="seed")
     network.size <- pmax(network.size,nr+!is.seed)
     data.title <- "reported network size\n (augmented by the number of recruits and the recruiter, if any)"
    }else{
     data.title <- "reported network sizes"
    }
   }else{
     network.size <- as.numeric(x$data)
   }
   
   remvalues <- is.na(network.size)
   if(!is.null(x$rectime)){
     network.size[!remvalues] <- (network.size[!remvalues])[order(x$rectime)]
   }
   rescale <- exp(ifelse(is.null(x$mem.optimism.prior),0,x$mem.optimism.prior))
   Kmax <- max(c(length(x$predictive.visibility)*rescale,network.size[!remvalues]))
   ns.prob <- tabulate(round(network.size[!remvalues]),nbins=Kmax) #, nbins=max(x$data))
   ns.prob <- ns.prob/sum(ns.prob)
#  med.vis <- which.max(cumsum(x$predictive.visibility)>=0.5)
   if(control$type %in% c("visibility","all") & x$visibility){
    aaa <- graphics::barplot(ns.prob,names.arg=1:Kmax,add=FALSE,axes=TRUE,width=rep(0.5,length(ns.prob)),space=1,col=0,
     xlab="visibility",ylab="probability", xlim=c(1,Kmax), ylim=c(0,max(c(ns.prob,x$predictive.visibility))),
     main="Posterior with sample visibility\n histogram overlaid (median matched)", cex.main=cex.main,
     sub="(posterior in red)")
#   graphics::lines(x=-0.25+median(network.size[!remvalues])*seq_along(x$predictive.visibility)/med.vis,y=x$predictive.visibility, type='h', col='red', lwd=2)
    if(!is.null(x$mem.optimism.prior)){
      graphics::lines(x=-0.25+exp(x$mem.optimism.prior)*seq_along(x$predictive.visibility),y=x$predictive.visibility, type='h', col='red', lwd=2)
    }else{
      graphics::lines(x=-0.25+exp(x$mem.optimism.prior1)*seq_along(x$predictive.visibility),y=x$predictive.visibility, type='h', col='red', lwd=2)
     }
   }
  
   if(control$type %in% c("degree","all")){
    if(methods::is(x$data,"rds.data.frame") & !is.null(x$visibilities)){
#    nplot <- min(nrow(x$vsample),200)
     nplot <- ceiling(max(2000 / sum(!is.na(network.size)), 1))
     dat <- data.frame(x=rep(network.size[!is.na(network.size)],rep(nplot,sum(!is.na(network.size)))),
                       y=as.vector(x$vsample[1:nplot,1:sum(!is.na(network.size))]))
     gfit <- scam::scam(y ~ s(x, bs="mpi"), family=poisson(link=log), data=dat)
     pfit <- predict(gfit, newdata=list(x=1:max(network.size,na.rm=TRUE)), type="response", se.fit=FALSE)
     graphics::plot(y=x$visibilities, x=network.size,ylab="estimated visibilities",
       xlab=data.title, ylim=range(c(x$visibilities,pfit),na.rm=TRUE),
       main="Estimated Visibilites for each respondent", cex.main=cex.main,log=log.degree)
     lines(x=1:max(network.size,na.rm=TRUE),y=pfit)
     qs <- apply(x$vsample,2,stats::quantile,probs=c(0.25,0.75))
     errorbar(y=x$visibilities, x=network.size,yminus=qs[1,],yplus=qs[2,],add=TRUE)
    }
   }
  }
 }
}else{
if(control$type %in% c("prior","all")){
  cy <- cumsum(lpriorm)
  #
  if(is.null(control$xlim)){control$xlim <- xp[which.max(cy>0.99)]}
  if(is.null(control$ylim)){control$ylim <- c(0,max(lpriorm))}
  graphics::plot(x=xp,y=lpriorm,type='l', xlab="population size", 
    main="prior for population size", cex.main=cex.main,
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
