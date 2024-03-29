% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.size.R
\name{summary.sspse}
\alias{summary.sspse}
\title{Summarizing Population Size Estimation Model Fits}
\usage{
\method{summary}{sspse}(object, support = 1000, HPD.level = 0.95, method = "bgk", ...)
}
\arguments{
\item{object}{an object of class \code{"sspse"}, usually, a result of a call
to \code{\link{posteriorsize}}.}

\item{support}{the number of equally-spaced points to use for the support of
the estimated posterior density function.}

\item{HPD.level}{numeric; probability level of the highest probability
density interval determined from the estimated posterior.}

\item{method}{character; The method to use for density estimation (default Gaussian Kernel; "bgk").
"Bayes" uses a Bayesian density estimator which has good properties.}

\item{\dots}{further arguments passed to or from other methods.}
}
\value{
The function \code{summary.sspse} computes and returns a two row matrix of
summary statistics of the prior and estimated posterior distributions. The rows correspond to the \code{Prior} and the
\code{Posterior}, respectively. 
The rows names are \code{Mean}, \code{Median}, \code{Mode}, \code{25\%}, \code{75\%}, and \code{90\%}.
These correspond to the distributional mean, median, mode, lower quartile, upper quartile and 90\% quantile, respectively.
}
\description{
This is the \code{summary} method for class \code{"sspse"} objects.
These objects encapsulate an estimate of the posterior distribution of
the population size based on data collected by Respondent Driven Sampling.
The approach approximates the RDS via the Sequential Sampling model of Gile
(2008).  As such, it is referred to as the Sequential Sampling - Population Size Estimate (SS-PSE).
It uses the order of selection of the sample to provide information
on the distribution of network sizes over the population members.
\code{summary} method for class \code{"sspse"}. posterior distribution of
the population size based on data collected by Respondent Driven Sampling.
The approach approximates the RDS via the Sequential Sampling model of Gile
(2008). As such, it is referred to as the Sequential Sampling - Population Size Estimate (SS-PSE).
It uses the order of selection of the sample to provide information
on the distribution of network sizes over the population members.
}
\details{
\code{print.summary.sspse} tries to be smart about formatting the
coefficients, standard errors, etc. and additionally gives
\sQuote{significance stars} if \code{signif.stars} is \code{TRUE}.

Aliased coefficients are omitted in the returned object but restored by the
\code{print} method.

Correlations are printed to two decimal places (or symbolically): to see the
actual correlations print \code{summary(object)$correlation} directly.
}
\examples{

data(fauxmadrona)
# Here interval=1 so that it will run faster. It should be higher in a 
# real application.
fit <- posteriorsize(fauxmadrona, median.prior.size=1000,
                                 burnin=20, interval=1, samplesize=100)
summary(fit)

}
\seealso{
The model fitting function \code{\link{posteriorsize}},
\code{\link{summary}}.
}
\keyword{models}
