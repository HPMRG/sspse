#' Summarizing Population Size Estimation Model Fits
#' 
#' This is the \code{summary} method for class \code{"sspse"} objects.
#' These objects encapsulate an estimate of the posterior distribution of
#' the population size based on data collected by Respondent Driven Sampling.
#' The approach approximates the RDS via the Sequential Sampling model of Gile
#' (2008).  As such, it is referred to as the Sequential Sampling - Population Size Estimate (SS-PSE).
#' It uses the order of selection of the sample to provide information
#' on the distribution of network sizes over the population members.
#' \code{summary} method for class \code{"sspse"}. posterior distribution of
#' the population size based on data collected by Respondent Driven Sampling.
#' The approach approximates the RDS via the Sequential Sampling model of Gile
#' (2008). As such, it is referred to as the Sequential Sampling - Population Size Estimate (SS-PSE).
#' It uses the order of selection of the sample to provide information
#' on the distribution of network sizes over the population members.
#' 
#' \code{print.summary.sspse} tries to be smart about formatting the
#' coefficients, standard errors, etc. and additionally gives
#' \sQuote{significance stars} if \code{signif.stars} is \code{TRUE}.
#' 
#' Aliased coefficients are omitted in the returned object but restored by the
#' \code{print} method.
#' 
#' Correlations are printed to two decimal places (or symbolically): to see the
#' actual correlations print \code{summary(object)$correlation} directly.
#' 
#' @aliases summary.sspse
#' @param object an object of class \code{"sspse"}, usually, a result of a call
#' to \code{\link{posteriorsize}}.
#' @param support the number of equally-spaced points to use for the support of
#' the estimated posterior density function.
#' @param HPD.level numeric; probability level of the highest probability
#' density interval determined from the estimated posterior.
#' @param method character; The method to use for density estimation (default Gaussian Kernel; "bgk").
#' "Bayes" uses a Bayesian density estimator which has good properties.
#' @param \dots further arguments passed to or from other methods.
#' @return The function \code{summary.sspse} computes and returns a two row matrix of
#' summary statistics of the prior and estimated posterior distributions. The rows correspond to the \code{Prior} and the
#' \code{Posterior}, respectively. 
#' The rows names are \code{Mean}, \code{Median}, \code{Mode}, \code{25\%}, \code{75\%}, and \code{90\%}.
#' These correspond to the distributional mean, median, mode, lower quartile, upper quartile and 90\% quantile, respectively. 
#' @seealso The model fitting function \code{\link{posteriorsize}},
#' \code{\link{summary}}.
#' 
#' @keywords models
#' @examples
#' 
#' data(fauxmadrona)
#' # Here interval=1 so that it will run faster. It should be higher in a 
#' # real application.
#' fit <- posteriorsize(fauxmadrona, median.prior.size=1000,
#'                                  burnin=20, interval=1, samplesize=100)
#' summary(fit)
#' 
#' @method summary sspse
#' @export
summary.sspse <- function(object, support=1000, HPD.level=0.95, method="bgk",...){
#summary.sspse <- function(object, ...){
  p.args <- as.list( sys.call() )[-c(1,2)]
  formal.args<-formals(sys.function())[-1]
# control <- list(support=1000,HPD.level=0.95)

  control<-list(samples=4000, burnin=1000)
  names.formal.args <- names(formal.args)
  names.formal.args <- names.formal.args[-match("...",names.formal.args)]
  for(arg in names.formal.args){ control[arg]<-list(get(arg)) }
  for(arg in names(p.args)){ control[arg]<-list(get(arg)) }

#suppressMessages(require(locfit, quietly=TRUE))
out <- object$sample
if(is.null(out)){
  object$n <- min(object$x)
  object$lpriorm <- log(object$lprior)
}
if(!is.null(out)){
  outN <- out[,"N"]
  #a=locfit( ~ lp(outN, nn=0.35, h=0, maxk=500))
  xp <- seq(object$n,object$maxN, length=control$support)
# a=locfit::locfit( ~ lp(outN,nn=0.5))
# posdensN <- predict(a, newdata=xp)
  if(method=="bgk" || !requireNamespace("densEstBayes", quietly = TRUE)){
    a=bgk_kde(outN,n=2^(ceiling(log((object$maxN-object$n))/log(2))),MIN=object$n,MAX=object$maxN)
    posdensN <- stats::spline(x=a[1,],y=a[2,],xout=xp)$y
  }else{
    a=densEstBayes::densEstBayes(outN,method="NUTS",
      control=densEstBayes::densEstBayes.control(range.x=c(object$n*0.95,object$maxN*1.05),#numBins=min(401,round(length(outN)/2)),
                                                 nKept=control$samples,nWarm=control$burnin))
    xTrang <- seq(-0.05, 1.05, length = length(xp))
    Xg <- cbind(1,xTrang)
    Zg <- .ZOSull(xTrang,intKnots=a$intKnots,range.x=c(-0.05,1.05))
    Cg <- cbind(Xg,Zg)
    betauMCMC <- a$stochaFitObj$betauMCMC
    etaHatMCMC <- crossprod(t(Cg),betauMCMC)
    posdensN <- exp(apply(etaHatMCMC, 1, mean))
  }
  posdensN <- control$support*posdensN / ((object$maxN-object$n)*sum(posdensN))
  # Next from coda
  # hpd <- HPDinterval(object$sample[,"N"])[1:2]
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
