\name{posteriorsize}
\alias{posteriorsize}
\title{Bayesian inference for the population size based on RDS data}
\description{
  \code{\link{posteriorsize}} is used to compute the marginal posterior distribution of the population size based on data collected by Respondent Driven Sampling. The model assumes an SPPS model for RDS. It uses order of selection of the sample to provide information on the distribution of degrees over the population 
members. The prior for the class sizes is multinomial given the population size. The probability that a randomly chosen person has degree \eqn{k} is modeled as a Waring distribution with specified (hyper)parameters. 
}
\usage{
posteriorsize(s, maxN = 4 * length(s), K = 2 * max(s), nk = tabulate(s, 
    nbin = K), n = length(s), N = 0.5 * maxN, mean.prior.degree = 7, 
    sd.prior.degree = 2.2, df.mean.prior = 1, df.sd.prior = 5, 
    muproposal = 0.01, sigmaproposal = 0.1, samplesize = 1000, 
    burnin = 100, interval = 1, burnintheta = 500, parallel = 1, 
    seed = NULL, verbose = TRUE)
}
\arguments{
  \item{s}{vector of integers; the vector of degrees from the RDS in order 
    they are recorded.}
  \item{maxN}{integer; maximum possible population size.  The prior
    for population size is uniform from 1 to MaxN}
  \item{K}{count; the maximum degree for an individual.}
  \item{nk}{vector; the vector of counts for the number of people in the sample
with degree k. This is usually computed from \eqn{s} automatically and not
usually specified by the user.}
  \item{n}{integer; sample size, that is, the length of \eqn{s},
  automatically computed and not usually specified by the user.}
  \item{N}{integer; initial guess for population sizesample size, that
  is, starting value for MCMC sampler. But it doesn't really matter.}
  \item{mean.prior.degree}{scalar; The mean degree for the prior
    distribution for a randomly chosen person. Prior is a
    Poisson-lognormal with this mean.}
  \item{sd.prior.degree}{scalar; The standard deviation of the degree for 
    the prior distribution. Prior is a
    Poisson-lognormal with this standard deviation.}
  \item{mu}{scalar; The mean for the prior distribution for a randomly chosen
    person's degree. The prior is a Waring distribution with this mean.
    By default for the mean is 5.}
  \item{rho}{scalar; The scaling parameter for the prior distribution for a randomly chosen
person's degree. The prior is a Waring distribution with this power-law scale.
    By default for the scale is 3 (i.e., just not a variance).}
  \item{M}{count; the number of Monte-Carlo samples of class sizes to draw to compute the marginal posterior. The default is 100000.}
  \item{parallel}{count; the number of parallel processes to run for the Monte-Carlo sample.  This uses PVM or MPI and requires the \code{snow} package and its dependencies to be installed.
The default is 1, that is not to use parallel processing.
  }
  \item{seed}{integer; random number integer seed.  Defaults to
   \code{NULL} to use whatever the state of the random number
   generator is at the time of the call.
  }
  \item{return.all}{logical; if this is
    \code{TRUE}, the function will return a long-list of supplementary
information in addition to the marginal posterior value. By default it only
returns the (log) marginal posterior value.}
  \item{verbose}{logical; if this is
    \code{TRUE}, the program will print out additional
    information, including goodness of fit statistics.
  }
}

\value{\code{\link{posteriorsize}} returns a list
  consisting of the following elements:
  \item{N}{vector; the vector of population sizes at which the marginal posterior (over the degree class sizes) has been computed.} 
  \item{logprob}{vector; the logarithm of the
    posterior probabiliy for each population size \eqn{N} marginalized over the number of people in the population in each of the degree classes.}
  \item{s}{vector; the vector of degrees from the RDS in order they are recorded.}
  \item{K}{count; the maximum degree for an individual.} 
  \item{mu}{scalar; The mean for the prior distribution for a randomly chosen
    person's degree. The prior is a Waring distribution with this mean.}
  \item{rho}{scalar; The scaling parameter for the prior distribution for a 
randomly chosen person's degree. The prior is a Waring distribution with this 
power-law scale.}
}

\references{

Gile, Krista J. and Handcock, Mark S. (2009)
\emph{Estimation of the size of a population based on Respondent-driven sampling} 
    Working Paper,
Nuffield College.

Handcock MS (2003{{b}}).
 {\pkg{degreenet}: Models for Skewed Count Distributions Relevant
  to Networks}.
 Statnet Project, Seattle, WA.
 Version 1.0, \url{http://statnet.org}.

Gile, Krista J. and Handcock, Mark S. (2009)
 {\pkg{size}: {A} Package to Estimate Population Size from Respondent Driven
Sampling Data}.
 Statnet Project, Seattle, WA.
 Version 1, \url{http://statnet.org}.

Handcock MS, Hunter DR, Butts CT, Goodreau SM, Morris M (2003{{b}}).
 {\pkg{statnet}: Software Tools for the Statistical Modeling of
  Network Data}.
 Statnet Project, Seattle, WA.
 Version 2, \url{http://statnet.org}.
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
# Create an RD  sample
#
set.seed(1)
pop<-sample(1:K, size=N0, replace = TRUE, prob = probs)
s <- rep(0,n)
rpop<-tabulate(pop, nbin=K)
for(i in 1:n){
 s[i] <- sample(1:K, size=1,prob=rpop*(1:K))
 rpop[s[i]] <- rpop[s[i]] - 1 
}
 
out <- posteriorsize(s=s,K=K,mu=mu,rho=rho,N=seq(180,320,length=10))
out
prob <- exp(out$logprob)/sum(exp(out$logprob))
pdf("bbb.pdf")
plot(x=out$N, y=prob)
lines(x=out$N, y=prob)
abline(v=N0)
abline(v=n,lty=2)
}
\keyword{models}
