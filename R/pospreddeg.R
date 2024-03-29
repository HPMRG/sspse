#' Compute the posterior predictive p-values for the reported network sizes
#' 
#' This function extracts from an
#' estimate of the posterior distribution of the
#' population size based on data collected by Respondent Driven Sampling. The
#' approach approximates the RDS via the Sequential Sampling model of Gile
#' (2008).  As such, it is referred to as the Sequential Sampling - Population Size Estimate (SS-PSE).
#' It uses the order of selection of the sample to provide information
#' on the distribution of network sizes over the population members.
#' 
#' It computes the posterior predictive distribution for each reported network size 
#' and computes the percentile rank of the reported network size within that posterior.
#' The percentile rank should be about 0.5 for a well specified model, but could be close to
#' uniform if there is little information about the reported network size.
#' The percentile ranks should not be extreme (e.g., close to zero or one) on a consistent basis
#' as this indicates a misspecified model.
#' 
#' @param x an object of class \code{"sspse"}, usually, a result of a call
#' to \code{oosteriorsize}.
#' @param order.by.recruitment.time logical; If \code{TRUE}, 
#' the reorder the input data by the recruitment time 
#' @seealso The model fitting function \code{\link{posteriorsize}},
#' \code{\link[graphics]{plot}}.
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
#' R package, Los Angeles, CA.  Version 0.5, \url{https://hpmrg.org/sspse/}.
#' 
#' Handcock MS (2003).  \pkg{degreenet}: Models for Skewed Count Distributions
#' Relevant to Networks.  Statnet Project, Seattle, WA.  Version 1.2,
#' \url{https://statnet.org/}.
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
#'\dontrun{
#' data(fauxmadrona)
#' # Here interval=1 so that it will run faster. It should be higher in a 
#' # real application.
#' fit <- posteriorsize(fauxmadrona, median.prior.size=1000,
#'                                  burnin=20, interval=1, samplesize=100)
#' summary(fit)
#' # Let's look at some MCMC diagnostics
#' pospreddeg(fit)
#' }
#' 
#' @export
pospreddeg <- function(x, order.by.recruitment.time=FALSE){
  deg <- RDS::get.net.size(x$data)
  if(order.by.recruitment.time){
    recruitment.time <- RDS::get.recruitment.time(x$data,  wave.fallback = TRUE)
    deg <- deg[order(recruitment.time)]
  }
  mem.optimism <- exp(x$mem.optimism.prior)
  mem.scale <- x$sample[,"mem.scale"]
  d <- 0:(10*x$K)
  ppd <- function(vs, d, deg, dlg){
   p <- matrix(0,ncol=ncol(vs), nrow=nrow(vs))
   for(i in 1:nrow(vs)){
    for(j in 1:ncol(vs)){
      ldraw <- d*log(vs[i,j])-vs[i,j]*mem.optimism-abs(d-mem.optimism*vs[i,j])/mem.scale[i] - dlg
      draw <- exp(ldraw) / sum(exp(ldraw))
      p[i,j] <- sum(draw[deg[j] > d])+0.5*draw[deg[j] == d]
    }
   }
   apply(p,2,mean)
  }
  pvalues <- ppd(x$vsample, d, deg, lgamma(d+1))
  class(pvalues) <- "pospreddeg"
  pvalues
}

#' Plots the posterior predictive p-values for the reported network sizes
#' 
#' This function extracts from an
#' estimate of the posterior distribution of the
#' population size based on data collected by Respondent Driven Sampling. The
#' approach approximates the RDS via the Sequential Sampling model of Gile
#' (2008).  As such, it is referred to as the Sequential Sampling - Population Size Estimate (SS-PSE).
#' It uses the order of selection of the sample to provide information
#' on the distribution of network sizes over the population members.
#' 
#' It computes the posterior predictive distribution for each reported network size 
#' and computes the percentile rank of the reported network size within that posterior.
#' The percentile rank should be about 0.5 for a well specified model, but could be close to
#' uniform if there is little information about the reported network size.
#' The percentile ranks should not be extreme (e.g., close to zero or one) on a consistent basis
#' as this indicates a misspecified model.
#' 
#' @param x an object of class \code{"pospreddeg"}, usually, a result of a call
#' to \code{pospreddeg}.
#' @param main character; title for the plot
#' @param nclass count; The number of classes for the histogram plot
#' @param hist logical; If \code{TRUE} plot a histogram of the p-values rather than
#' a density estimate.
#' @param ylim two-vector; lower and upper limits of vertical/density axis.
#' @param order.by.recruitment.time logical; If \code{TRUE}, 
#' the reorder the input data by the recruitment time 
#' @param \dots further arguments passed to or from other methods.
#' @seealso The model fitting function \code{\link{posteriorsize}},
#' \code{\link[graphics]{plot}}.
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
#' R package, Los Angeles, CA.  Version 0.5, \url{https://hpmrg.org/sspse/}.
#' 
#' Handcock MS (2003).  \pkg{degreenet}: Models for Skewed Count Distributions
#' Relevant to Networks.  Statnet Project, Seattle, WA.  Version 1.2,
#' \url{https://statnet.org/}.
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
#'\dontrun{
#' data(fauxmadrona)
#' # Here interval=1 so that it will run faster. It should be higher in a 
#' # real application.
#' fit <- posteriorsize(fauxmadrona, median.prior.size=1000,
#'                                  burnin=10, interval=1, samplesize=50)
#' summary(fit)
#' # Let's look at some MCMC diagnostics
#' plot(pospreddeg(fit))
#' }
#' 
#' @method plot pospreddeg
#' @export
plot.pospreddeg <- function(x,
    main="Posterior Predictive p-values for the self-reported network sizes", 
    nclass=20, hist=FALSE, ylim=c(0,2), order.by.recruitment.time=FALSE, ...){
  if(hist){
   outp <- hist(x, nclass=nclass, xlab="posterior p-value", 
     ylab="Density",xlim=c(0,1), ylim=ylim, probability=TRUE,
     main=main, ...)
   abline(h=1,lty=2)
  }else{
   xTrang <- seq(0, 1, length = 100)
   posdensp=KernSmooth::bkde(x=x, kernel = "normal", gridsize = length(xTrang),
            range.x=c(0,1))
   xTrang <- posdensp$x
   posdensp <- posdensp$y
   #
   outp <- graphics::plot(x=xTrang,y=posdensp,type='l', xlab="posterior p-value", 
     ylab="Density",xlim=c(0,1), ylim=ylim,
     main=main, ...)
   abline(h=1,lty=2)
  }
}
