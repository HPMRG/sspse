\name{margposN}
\alias{margposN}
\title{Inference for the population size based on RDS data}
\description{
  \code{\link{margposN}} is used to compute the marginal posterior distribution of the population size based on data collected by Respondent Driven Sampling. The model assumes an SPPS model for RDS. It uses order of selection of the sample to provide information on the distribution of degrees over the population 
members. The prior for the class sizes is multinomial given the population size. The probability that a randomly chosen person has degree \eqn{k} is modeled as a Waring distribution with specified (hyper)parameters. 
}
\usage{
margposN(N,K,s,mu=5,rho=3,M=100000,
         parallel=1,seed=NULL,verbose=FALSE, 
         n=tabulate(s,nbin=K),
         return.all=FALSE)
}
\arguments{
  \item{N}{count; the population size to compute the marginal posterior (over
the degree class sizes).} 
  \item{K}{count; the maximum degree for an individual.} 
  \item{s}{vector; the vector of degrees from the RDS in order they are
recorded.}
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

\value{\code{\link{margposN}} returns a scalar which is the logarithm of the
posterior of the population size \eqn{N} marginalized over the number of people
in the population in each of the degree classes.
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

N <- 500
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
pop<-sample(1:K, size=N, replace = TRUE, prob = probs)
s <- rep(0,n)
rpop<-tabulate(pop, nbin=K)
for(i in 1:n){
 s[i] <- sample(1:K, size=1,prob=rpop*(1:K))
 rpop[s[i]] <- rpop[s[i]] - 1 
}

margposN(N=N,K=K,s=s,mu=mu,rho=rho,M=10000)
}
\keyword{models}
