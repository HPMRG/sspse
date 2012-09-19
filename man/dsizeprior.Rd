\name{dsizeprior}
\alias{dsizeprior}
\title{Prior distributions for the size of a hidden population}
\description{
\code{\link{dsizeprior}} computes the 
prior distribution of the population size of a hidden population.
The prior is intended to be used in Bayesian inference for the population size
based on data
collected by Respondent Driven Sampling.
}
\usage{
 dsizeprior(n,
   type=c("proportion","nbinom","pln","flat","continuous"),
   mean.prior.size=NULL, sd.prior.size=NULL,
   mode.prior.sample.proportion=0.5,
   median.prior.sample.proportion=NULL,
   median.prior.size=NULL,
   median.mid.prior.size=NULL,
   mode.prior.size=NULL,
   quartiles.prior.size=NULL,
   effective.prior.df=1,
   alpha = NULL,
   maxN=NULL,
   log=FALSE,
   maxbeta = 100,
   maxNmax = 2e+05,
   verbose=TRUE)
}
\arguments{
  \item{n}{count; the sample size.}
  \item{type}{character; the type of parametric distribution to use for the prior on population size. The options are \code{proportion}
(for a prior on the sample proportion (i.e. \eqn{n/N}), 
\code{nbinom} (Negative-Binomial), \code{pln} (Poisson-log-normal),
\code{flat} (uniform), and \code{continuous} (the continuous version of
\code{proportion}).
The default is \code{proportion}.}
  \item{mean.prior.size}{scalar; A hyperparameter being the mean of the prior distribution on the population size.}
  \item{sd.prior.size}{scalar; A hyperparameter being the standard deviation of the prior distribution on the population size.}
  \item{mode.prior.sample.proportion}{scalar; A hyperparameter being the mode of the prior distribution on the sample proportion \eqn{n/N}.}
  \item{median.prior.sample.proportion}{scalar; A hyperparameter being the median of the prior distribution on the sample proportion \eqn{n/N}.}
  \item{median.prior.size}{scalar; A hyperparameter being the mode of the prior distribution on the population size.}
  \item{mode.prior.size}{scalar; A hyperparameter being the mode of the prior distribution on the population size.}
  \item{quartiles.prior.size}{vector of length 2; A pair of hyperparameters being the lower and upper quartiles of the prior distribution on the population size.
For example, \code{quartiles.prior.size=c(1000,4000)} corresponds to a prior
where the lower quartile (25\%) is 1000 and the upper (75\%) is 4000.}
  \item{effective.prior.df}{scalar; A hyperparameter being the effective number of samples worth of information represented in the prior distribution on the population size. By default this is 1, but it can be greater (or less!) to allow for different levels of uncertainty.}
  \item{alpha}{scalar; A hyperparameter being the first parameter of the beta prior model for the sample proportion. By default this is NULL, meaning that 1 is chosen. it can be any value at least 1 to allow for different levels of uncertainty.}
  \item{maxN}{integer; maximum possible population size. By default this is determined from an upper quantile of the prior distribution.}
  \item{log}{logical; return the prior or the the logarithm of the prior.}
  \item{maxbeta}{integer; maximum beta in the prior for population size. By default this is determined to ensure numerical stability.}
  \item{maxNmax}{integer; maximum possible population size. By default this is determined to ensure numerical stability.}
  \item{verbose}{logical; if this is
    \code{TRUE}, the program will print out additional
    information, including goodness of fit statistics.
  }
}

\value{\code{\link{dsizeprior}} returns a list
  consisting of the following elements:
  \item{x}{vector; vector of degrees \code{1:N} at which the prior PMF is
computed.}
  \item{lpriorm}{vector; vector of probabilities corresponding to the values in \code{x}.}
  \item{N}{scalar; a starting value for the population size computed from the
prior.}
  \item{maxN}{integer; maximum possible population size. By default this is determined from an upper quantile of the prior distribution.}
  \item{mean.prior.size}{scalar; A hyperparameter being the mean of the prior distribution on the population size.}
% \item{sd.prior.size}{scalar; A hyperparameter being the standard deviation of the prior distribution on the population size.}
  \item{mode.prior.size}{scalar; A hyperparameter being the mode of the prior distribution on the population size.}
  \item{effective.prior.df}{scalar; A hyperparameter being the effective number of samples worth of information represented in the prior distribution on the population size. By default this is 1, but it can be greater (or less!) to allow for different levels of uncertainty.}
  \item{mode.prior.sample.proportion}{scalar; A hyperparameter being the mode of the prior distribution on the sample proportion \eqn{n/N}.}
  \item{median.prior.size}{scalar; A hyperparameter being the mode of the prior distribution on the population size.}
  \item{beta}{scalar; A hyperparameter being the second parameter of the Beta distribuion that is a component of the prior distribution on the sample proportion \eqn{n/N}.}
  \item{type}{character; the type of parametric distribution to use for the prior on population size. The options are \code{proportion}
(for a prior on the sample proportion (i.e. \eqn{n/N}), 
\code{nbinom} (Negative-Binomial), \code{pln} (Poisson-log-normal),
\code{flat} (uniform), and \code{continuous} (the continuous version of
\code{proportion}).
The default is \code{proportion}.}
}

\section{Details on priors}{
The best way to specify the prior is via the hyperparameter \code{mode.prior.size}
which specifies the mode of the prior distribution on the population size.
You can alternatively specify the hyperparameter \code{median.prior.size}
which specifies the median of the prior distribution on the population size,
or \code{mode.prior.sample proportion} which specifies the mode of the prior
distribution on the proportion of the population size in the sample.
}

\references{

Gile, Krista J. (2008)
\emph{Inference from Partially-Observed Network Data},
Ph.D. Thesis, Department of Statistics, University of Washington.

Gile, Krista J. and Handcock, Mark S. (2010)
\emph{Respondent-Driven Sampling: An Assessment of Current Methodology},
To appear, Sociological methodology, 40. 

Gile, Krista J. and Handcock, Mark S. (2009)
 \pkg{size}: A Package to Estimate Population Size from Respondent Driven
Sampling Data.
 Statnet Project, Seattle, WA.
 Version 1, \url{http://statnet.org}.

Handcock MS (2003b).
 \pkg{degreenet}: Models for Skewed Count Distributions Relevant
  to Networks.
 Statnet Project, Seattle, WA.
 Version 1.0, \url{http://statnetproject.org}.

Handcock, Mark S.,  Gile, Krista J. and Mar, Corinne M. (2009)
\emph{Estimating Hidden Population Size using Respondent-Driven Sampling
Data}, manuscript.

Handcock MS, Hunter DR, Butts CT, Goodreau SM, Morris M (2003b).
 \pkg{statnet}: Software Tools for the Statistical Modeling of
  Network Data.
 Statnet Project, Seattle, WA.
 Version 2, \url{http://statnetproject.org}.

Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b).
 \pkg{ergm}: A Package to Fit, Simulate and Diagnose
  Exponential-Family Models for Networks.
 \emph{Journal of Statistical Software}, 24(3).
 \url{http://www.jstatsoft.org/v24/i03/}.
}

\seealso{network, statnet, degreenet}

\examples{
dsizeprior(n=100,
           type="proportion",
           mode.prior.size=1000)
}
\keyword{models}
