\name{posteriorsize}
\alias{posteriorsize}
%\title{Bayesian inference for the population size based on RDS data}
\title{Estimating hidden population size using RDS data}
\description{
\code{\link{posteriorsize}} computes the 
posterior distribution of the population size based on data
collected by Respondent Driven Sampling. The approach approximates the
RDS via the Sequential Sampling model of Gile (2008). 
It uses the order of selection of the sample to
provide information on the distribution of network sizes 
over the population members. 

It does not use the disease states of the sampled individuals to improve the
estimation.  See \code{\link{posteriordisease}} for a version that uses the
sample of degrees only. 

%By default, the probability that a randomly chosen
%person has degree \eqn{k} is modeled as a Conway-Maxwell-Poisson distribution
%with specified (hyper)parameters.
}
\usage{
 posteriorsize(s,
                  mean.prior.degree=7,
                  sd.prior.degree=3,
                  df.mean.prior=1,df.sd.prior=5,
                  Np=0,
                  samplesize=1000,burnin=100,interval=1,burnintheta=500,
                  priorsizedistribution=c("proportion","nbinom","pln","flat","continuous"),
                  mean.prior.size=NULL, sd.prior.size=NULL,
                  mode.prior.sample.proportion=0.5,
                  median.prior.size=NULL,
                  mode.prior.size=NULL,
                  degreedistribution=c("cmp","nbinom","pln"),
                  maxN=NULL,
                  K=2*max(s), n=length(s),
                  nk=tabulate(s,nbin=K),
                  muproposal=0.1,
                  sigmaproposal=0.15,
                  parallel=1, seed=NULL,
                  verbose=TRUE)
}
\arguments{
  \item{s}{vector of integers; the vector of degrees from the RDS in order 
    they are recorded.}
  \item{mean.prior.degree}{scalar; A hyper parameter being the mean degree for the prior
    distribution for a randomly chosen person. The prior has this
    mean.}
  \item{sd.prior.degree}{scalar;  A hyper parameter being the standard deviation of the degree for a randomly chosen person.  The prior has this standard deviation.}
  \item{df.mean.prior}{scalar;  A hyper parameter being the degrees-of-freedom of the prior for the mean. This gives the equivalent sample size that would contain the same amount of information inherent in the prior.}
  \item{df.sd.prior}{scalar;  A hyper parameter being the degrees-of-freedom of the prior for the standard deviation. This gives the equivalent sample size that would contain the same amount of information inherent in the prior for the standard deviation.}
  \item{Np}{integer; The overall degree distribution 
  is a mixture of the \code{1:Np} rates and a parametric
  degree distribution model truncated below Np. Thus the model
  fits the proportions of the
population with degree \code{1:Np} each with a separate parameter. This
should adjust for an lack-of-fit of the parametric degree
distribution model  at lower degrees, although it also changes
the model away from the parametric degree distribution model.}
  \item{samplesize}{count; the number of Monte-Carlo samples to draw to compute the posterior. This is the number returned by the Metropolis-Hastings algorithm.The default is 1000.}
  \item{burnin}{count; the number of proposals before any MCMC sampling
    is done. It typically is set to a fairly large number.}
  \item{interval}{count; the number of proposals between sampled statistics.}
  \item{burnintheta}{count; the number of proposals in the Metropolis-Hastings
sub-step for the degree distribution parameters (\eqn{\theta}) before any MCMC sampling is done. It typically is set to a modestly large number.}
  \item{priorsizedistribution}{character; the type of parametric distribution to use for the prior on population size. The options are \code{proportion}
(for a prior on the sample proportion (i.e. \eqn{n/N}), 
\code{nbinom} (Negative-Binomial), \code{pln} (Poisson-log-normal),
\code{flat} (uniform),
and \code{continuous} (the continuous version of \code{proportion}).The default is \code{proportion}.
The default is \code{proportion}.}
  \item{mean.prior.size}{scalar; A hyperparameter being the mean of the prior distribution on the population size.}
  \item{sd.prior.size}{scalar; A hyperparameter being the standard deviation of the prior distribution on the population size.}
  \item{mode.prior.sample.proportion}{scalar; A hyperparameter being the mode of the prior distribution on the sample proportion \eqn{n/N}.}
  \item{median.prior.size}{scalar; A hyperparameter being the mode of the prior distribution on the population size.}
  \item{mode.prior.size}{scalar; A hyperparameter being the mode of the prior distribution on the population size.}
  \item{degreedistribution}{count; the parametric distribution to
  use for the individual network sizes (i.e., degrees). The
  options are \code{cmp}, \code{nbinom}, and \code{pln}.
  These correspond to the Conway-Maxwell-Poisson,
  Negative-Binomial, and Poisson-log-normal. The default is \code{cmp}.}
  \item{maxN}{integer; maximum possible population size. By default this is determined from an upper quantile of the prior distribution.}
  \item{K}{count; the maximum degree for an individual. This is usually
calculated as twice the maximum observed degree.}
  \item{nk}{vector; the vector of counts for the number of people in the sample
with degree k. This is usually computed from \eqn{s} automatically and not
usually specified by the user.}
  \item{muproposal}{scalar; The standard deviation of the proposal distribution for the mean degree.}
  \item{sigmaproposal}{scalar; The standard deviation of the proposal
distribution for the standard deviation of the degree.}
  \item{parallel}{count; the number of parallel processes to run for the Monte-Carlo sample.  This uses PVM or MPI and requires the \code{snow} package and its dependencies to be installed.
The default is 1, that is not to use parallel processing.}
  \item{seed}{integer; random number integer seed.  Defaults to
   \code{NULL} to use whatever the state of the random number
   generator is at the time of the call.
  }
  \item{verbose}{logical; if this is
    \code{TRUE}, the program will print out additional
    information, including goodness of fit statistics.
  }
}

\value{\code{\link{posteriorsize}} returns a list
  consisting of the following elements:
  \item{pop}{vector; The final posterior draw for the degrees of the population. The first \eqn{n} at the sample and the reminder are non-sequenced. The length of the output vector is \code{maxN}, but only the first \code{N} degrees are relevant (the vector is padded by zeros.}
  \item{K}{count; the maximum degree for an individual. This is usually
calculated as twice the maximum observed degree.}
  \item{n}{count; the sample size.}
  \item{samplesize}{count; the number of Monte-Carlo samples to draw to compute the posterior. This is the number returned by the Metropolis-Hastings algorithm.The default is 1000.}
  \item{burnin}{count; the number of proposals before any MCMC sampling
    is done. It typically is set to a fairly large number.}
  \item{interval}{count; the number of proposals between sampled statistics.}
  \item{mu}{scalar; The hyper parameter \code{mean.prior.degree} 
being the mean degree for the prior distribution for a randomly chosen person. The prior has this mean.}
  \item{sigma}{scalar;  The hyper parameter \code{sd.prior.degree} being the standard deviation of the degree for 
    a randomly chosen person. The prior has this standard deviation.}
  \item{df.sd.prior}{scalar;  A hyper parameter being the degrees-of-freedom of the prior for the standard deviation. This gives the equivalent sample size that would contain the same amount of information inherent in the prior for the standard deviation.}
  \item{Np}{integer; The overall degree distribution
  is a mixture of the \code{1:Np} rates and a parametric
  degree distribution model truncated below Np. Thus the model
  fits the proportions of the
population with degree \code{1:Np} each with a separate parameter. This
should adjust for an lack-of-fit of the parametric degree
distribution model  at lower degrees, although it also changes
the model away from the parametric degree distribution model.}
  \item{muproposal}{scalar; The standard deviation of the proposal distribution for the mean degree.}
  \item{sigmaproposal}{scalar; The standard deviation of the proposal
distribution for the standard deviation of the degree.}
  \item{N}{vector of length 5; summary statistics for the posterior population size. 
   \describe{
      \item{\code{MAP}}{maximum aposteriori value of N}
      \item{\code{Mean AP}}{mean aposteriori value of N}
      \item{\code{P025}}{the 2.5th percentile of the (posterior) distribution
for the N. That is, the lower point on a 95\% probability interval.}
      \item{\code{P975}}{the 97.5th percentile of the (posterior) distribution
for the N. That is, the upper point on a 95\% probability interval.}
   }
  }
  \item{maxN}{integer; maximum possible population size. By default this is determined from an upper quantile of the prior distribution.}
 \item{sample}{matrix of dimension \code{samplesize}\eqn{\times} \code{10} matrix of summary statistics from the posterior. this is also an object of class \code{mcmc} so it can be plotted and summarized via the \code{mcmc.diagnostics} function in the \code{ergm} package (and also the \code{coda} package). The statistics are:
   \describe{
      \item{\code{N}}{population size.}
      \item{mu}{scalar; The mean degree for the prior distribution for a randomly chosen person. The prior has this mean.}
      \item{sigma}{scalar; The standard deviation of the degree for 
    a randomly chosen person. The prior has this standard deviation.}
      \item{degree1}{scalar; the number of nodes of degree 1 in the population (it is assumed all nodes have degree 1 or more).}
      \item{lambda}{scalar; This is only present for the \code{cmp} model. It is the \eqn{\lambda} parameter in the standard parametrization of the Conway-Maxwell-Poisson model for the degree distribution.}
      \item{nu}{scalar; This is only present for the \code{cmp} model. It is the \eqn{\nu} parameter in the standard parametrization of the Conway-Maxwell-Poisson model for the degree distribution.}
     }
}
  \item{lpriorm}{vector; the vector of (log) prior probabilities on each value
of \eqn{m=N-n} - that is, the number of unobserved members of the population.}
  \item{burnintheta}{count; the number of proposals in the Metropolis-Hastings
sub-step for the degree distribution parameters (\eqn{\theta}) before any MCMC sampling is done. It typically is set to a modestly large number.}
  \item{fVerbose}{logical; if this is
    \code{TRUE}, the program printed out additional
    information, including goodness of fit statistics.}
  \item{predictive.degree.count}{vector; a vector of length the maximum degree (\code{K}) (by default \code{K=2*max(sample degree)}). 
  The \code{k}th entry is the posterior predictive number persons with degree \code{k}.
    That is, it is the posterior predictive distribution of the number of
people with each degree in the population.}
  \item{predictive.degree}{vector; a vector of length the maximum degree (\code{K}) (by default \code{K=2*max(sample degree)}). 
  The \code{k}th entry is the posterior predictive proportion of persons with degree \code{k}.
    That is, it is the posterior predictive distribution of the proportion of
people with each degree in the population.}
 \item{MAP}{vector of length 10 of MAP estimates corresponding to the output \code{sample}. These are:
   \describe{
      \item{\code{N}}{population size.}
      \item{mu}{scalar; The mean degree for the prior distribution for a randomly chosen person. The prior has this mean.}
      \item{sigma}{scalar; The standard deviation of the degree for 
    a randomly chosen person. The prior has this standard deviation.}
      \item{degree1}{scalar; the number of nodes of degree 1 in the population (it is assumed all nodes have degree 1 or more).}
      \item{lambda}{scalar; This is only present for the \code{cmp} model. It is the \eqn{\lambda} parameter in the standard parametrization of the Conway-Maxwell-Poisson model for the degree distribution.}
      \item{nu}{scalar; This is only present for the \code{cmp} model. It is the \eqn{\nu} parameter in the standard parametrization of the Conway-Maxwell-Poisson model for the degree distribution.}
     }
  }
  \item{mode.prior.sample.proportion}{scalar; A hyperparameter being the mode of the prior distribution on the sample proportion \eqn{n/N}.}
  \item{median.prior.size}{scalar; A hyperparameter being the mode of the prior distribution on the population size.}
  \item{mode.prior.size}{scalar; A hyperparameter being the mode of the prior distribution on the population size.}
% \item{mean.prior.size}{scalar; A hyperparameter being the mean of the prior distribution on the population size.}
  \item{degreedistribution}{count; the parametric distribution to
  use for the individual network sizes (i.e., degrees). The
  options are \code{cmp}, \code{nbinom}, and \code{pln}.
  These correspond to the Conway-Maxwell-Poisson,
  Negative-Binomial, and Poisson-log-normal. The default is \code{cmp}.}
  \item{priorsizedistribution}{character; the type of parametric distribution to use for the prior on population size. The options are \code{proportion}
(for a prior on the sample proportion (i.e. \eqn{n/N}), 
\code{nbinom} (Negative-Binomial), \code{pln} (Poisson-log-normal),
\code{flat} (uniform),
and \code{continuous} (the continuous version of \code{proportion}).The default is \code{proportion}.
The default is \code{proportion}.
}

% \item{\dots}{Additional
%   arguments, to be passed to lower-level functions in the future.
% }
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
 {\pkg{size}: {A} Package to Estimate Population Size from Respondent Driven
Sampling Data}.
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

\seealso{network, statnet, snow, degreenet}

\examples{
library(degreenet)

N0 <- 200
n <- 100
K <- 10
rho <- 3
mu <- 5

# Look at the degree distribution for the prior
p0 <- mu*(rho-2)-rho+1
probs <- dwar(v=c(rho,p0),x=1:K)
probs <- probs / sum(probs)
# Plot these if you want
# plot(x=1:K,y=probs,type="l")
# points(x=1:K,y=probs)
#
# Create an RD sample
#
set.seed(1)
pop<-sample(1:K, size=N0, replace = TRUE, prob = probs)
s<-sample(pop, size=n, replace = FALSE, prob = pop)
dis<-sample(1*(pop>quantile(pop,0.2)), size=n, replace = TRUE, prob = pop)
 
out <- posteriorsize(s=s,dis=dis)
out
prob <- exp(out$logprob)/sum(exp(out$logprob))
x <- n+seq(along=out$lpriorm)
#pdf("bbb.pdf")
plot(density(out$sample[,"N"]),
     xlab="population size",ylab="Probability", type="n")
lines(x=x, y=exp(out$lpriorm), lty=2)
abline(v=N0)
abline(v=n,lty=2)
}
\keyword{models}
