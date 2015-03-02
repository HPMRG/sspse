\name{plot.psess}
\alias{plot.psess}
\title{Plot Summary and Diagnostics for Population Size Estimation Model Fits}
\usage{
\method{plot}{psess}(x,xlim=NULL,data=NULL,support=1000,
        HPD.level=0.90,N=NULL,ylim=NULL,mcmc=FALSE,type="both", \ldots)
}
\arguments{
  \item{x}{an object of class \code{"plot.psess"}, usually, a result of a
    call to \code{plot.psess}.}
  \item{xlim}{the (optional) x limits (x1, x2) of the plot of the posterior of the population size.}
  \item{data}{Optionally, the vector of degrees from the RDS in order they are recorded and as
passed to \code{\link{posteriorsize}}.}
  \item{support}{the number of equally-spaced points to use for the support of the estimated
        posterior density function.}
  \item{HPD.level}{numeric; probability level of the highest probability density interval
    determined from the estimated posterior.}
  \item{N}{Optionally, an estimate of the population size to mark on the plots as a reference
point.}
  \item{ylim}{the (optional) vertical limits (y1, y2) of the plot of the posterior of the population size. A vertical axis is the probability density scale.}
  \item{mcmc}{logical; If TRUE, additionally create simple diagnostic plots for the MCMC sampled statistics produced from the fit.}
  \item{type}{character; This controls the types of plots produced.
 If \code{"N"}, a density plot of the posterior for population size is produced.
 If \code{"others"}, a density plot of the prior for population size,
 a density plot of the posterior for mean network size in the population,
 the posterior for standard deviation of the network size, and a density plot of 
 the posterior mean network size distribution with sample histogram overlaid is produced.
 If \code{"both"}, then all plots for \code{"N"} and \code{"others"} are produced.
  }
  \item{\dots}{further arguments passed to or from other methods.}
}
\description{
\code{plot} method for class \code{"psess"}.
posterior distribution of the population size based on data
collected by Respondent Driven Sampling. The approach approximates the
RDS via the Sequential Sampling model of Gile (2008).
It uses the order of selection of the sample to
provide information on the distribution of network sizes
over the population members.

}
\details{
Six plots (selectable by \code{which}) are currently available: a plot
  of residuals against fitted values, a Scale-Location plot of
  \eqn{\sqrt{| residuals |}}{sqrt(| residuals |)}
  against fitted values, a Normal Q-Q plot, a
  plot of Cook's distances versus row labels, a plot of residuals
  against leverages, and a plot of Cook's distances against
  leverage/(1-leverage).  By default, the first three and \code{5} are provided.
}
\value{
  The function \code{plot.psess} computes and returns a list of plot
  statistics of the fitted linear model given in \code{object}, using
  the components (list elements) \code{"call"} and \code{"terms"}
  from its argument, plus
  \item{residuals}{the \emph{weighted} residuals, the usual residuals
    rescaled by the square root of the weights specified in the call to
    \code{psess}.}
  \item{coefficients}{a \eqn{p \times 4}{p x 4} matrix with columns for
    the estimated coefficient, its standard error, t-statistic and
    corresponding (two-sided) p-value.  Aliased coefficients are omitted.}
  \item{aliased}{named logical vector showing if the original
    coefficients are aliased.}
  \item{sigma}{the square root of the estimated variance of the random
    error
    \deqn{\hat\sigma^2 = \frac{1}{n-p}\sum_i{w_i R_i^2},}{\sigma^2 = 1/(n-p) Sum(w[i] R[i]^2),}
    where \eqn{R_i}{R[i]} is the \eqn{i}-th residual, \code{residuals[i]}.}
  \item{df}{degrees of freedom, a 3-vector \eqn{(p, n-p, p*)}, the first
    being the number of non-aliased coefficients, the last being the total
    number of coefficients.}
  \item{fstatistic}{(for models including non-intercept terms)
    a 3-vector with the value of the F-statistic with
    its numerator and denominator degrees of freedom.}
  \item{r.squared}{\eqn{R^2}, the \sQuote{fraction of variance explained by
    the model},
    \deqn{R^2 = 1 - \frac{\sum_i{R_i^2}}{\sum_i(y_i- y^*)^2},}{R^2 = 1 - Sum(R[i]^2) / Sum((y[i]- y*)^2),}
    where \eqn{y^*}{y*} is the mean of \eqn{y_i}{y[i]} if there is an
    intercept and zero otherwise.}
  \item{adj.r.squared}{the above \eqn{R^2} statistic
    \sQuote{\emph{adjusted}}, penalizing for higher \eqn{p}.}
  \item{cov.unscaled}{a \eqn{p \times p}{p x p} matrix of (unscaled)
    covariances of the \eqn{\hat\beta_j}{coef[j]}, \eqn{j=1, \dots, p}.}
  \item{correlation}{the correlation matrix corresponding to the above
    \code{cov.unscaled}, if \code{correlation = TRUE} is specified.}
  \item{symbolic.cor}{(only if \code{correlation} is true.)  The value
    of the argument \code{symbolic.cor}.}
  \item{na.action}{from \code{object}, if present there.}
}
\seealso{
  The model fitting function \code{\link{posteriorsize}}, \code{\link{plot}}.

  Function \code{\link{coef}} will extract the matrix of coefficients
  with standard errors, t-statistics and p-values.
}
\examples{
\donttest{
if(require("degreenet")){
N0 <- 200
n <- 180
K <- 10
rho <- 3
mu <- 5

p0 <- mu*(rho-2)-rho+1
probs <- dwar(v=c(rho,p0),x=1:K)
probs <- probs / sum(probs)
sum(probs*(1:K))
sqrt(sum(probs*(1:K)^2)-sum(probs*(1:K))^2)
pop<-sample(1:K, size=N0, replace = TRUE, prob = probs)
# The population is here is synthetic, being N=200
# draws from a Waring distribution with mena 5 and
# scaling parameter rho=3. The distribution is truncated at 10.

# Let's draw a sample
s<-sample(1:N0, size=n, replace = FALSE, prob = pop)
  
out <- posteriorsize(s=pop[s],K=K+5,interval=100,burnin=1000,
	 samplesize=100,
         priorsizedistribution="flat",
         maxN=500,maxbeta=30,verbose=TRUE)
# Let's look at some MCMC diagnostics
plot(out, HPD.level=0.9,mcmc=TRUE)
plot(out, HPD.level=0.9,data=pop[s])
summary(out, HPD.level=0.9)
}
}
}
\keyword{hplot}
