\name{margpossize}
\alias{margpossize}
\title{Inference for the population size based on RDS data}
\description{
  \code{\link{margpossize}} is used to compute the marginal posterior distribution of the population size based on data collected by Respondent Driven Sampling. The model assumes an SPPS model for RDS. It uses order of selection of the sample to provide information on the distribution of degrees over the population 
members. The prior for the class sizes is multinomial given the population size. The probability that a randomly chosen person has degree \eqn{k} is modeled as a Waring distribution with specified (hyper)parameters. 
}
\usage{
margpossize(s,N=trunc(length(s)*seq(1.1,4,length=10)+1),K=(1:max(s)),
          mu=5,rho=3,M=100000,
          parallel=1,seed=NULL, verbose=FALSE,
          n=tabulate(s,nbin=K),
          return.all=FALSE)
}
\arguments{
  \item{s}{vector; the vector of degrees from the RDS in order they are recorded.}
  \item{N}{vector; a vector of population sizes to compute the marginal posterior (over the degree class sizes) for.} 
  \item{K}{count; the maximum degree for an individual.} 
  \item{mu}{scalar; The mean for the prior distribution for a randomly chosen
    person's degree. The prior is a Waring distribution with this mean.
    By default for the mean is 5.}
  \item{rho}{scalar; The scaling parameter for the prior distribution for a randomly chosen
person's degree. The prior is a Waring distribution with this power-law scale.
    By default for the scale is 3 (i.e., just not a variance).}
  \item{M}{count; the number of Monte-Carlo samples of class sizes to draw to compute the marginal posterior. The default is 100000.
  }
  \item{parallel}{count; the number of parallel processes to run for the Monte-Carlo sample.  This uses PVM or MPI and requires the \code{snow} package and its dependencies to be installed.
The default is 1, that is not to use parallel processing.
  }
  \item{seed}{integer; random number integer seed.  Defaults to
   \code{NULL} to use whatever the state of the random number
   generator is at the time of the call.
  }
  \item{n}{vector; the vector of counts for the number of people in the sample
with degree k. This is usually computed from \eqn{s} automatically and not
usually specified by the user.}
  \item{return.all}{logical; if this is
    \code{TRUE}, the function will return a long-list of supplementary
information in addition to the marginal posterior value. By default it only
returns the (log) marginal posterior value.}
  \item{verbose}{logical; if this is
    \code{TRUE}, the program will print out additional
    information, including goodness of fit statistics.
  }
}

\value{\code{\link{margpossize}} returns a list
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

out <- margpossize(s=s,K=K,mu=mu,rho=rho,N=seq(180,320,length=10))
out
prob <- exp(out$logprob)/sum(exp(out$logprob))
pdf("bbb.pdf")
plot(x=out$N, y=prob)
lines(x=out$N, y=prob)
abline(v=N0)
abline(v=n,lty=2)
}
\keyword{models}
