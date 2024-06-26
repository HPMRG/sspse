% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pospreddeg.R
\name{pospreddeg}
\alias{pospreddeg}
\title{Compute the posterior predictive p-values for the reported network sizes}
\usage{
pospreddeg(x, order.by.recruitment.time = FALSE)
}
\arguments{
\item{x}{an object of class \code{"sspse"}, usually, a result of a call
to \code{oosteriorsize}.}

\item{order.by.recruitment.time}{logical; If \code{TRUE}, 
the reorder the input data by the recruitment time}
}
\description{
This function extracts from an
estimate of the posterior distribution of the
population size based on data collected by Respondent Driven Sampling. The
approach approximates the RDS via the Sequential Sampling model of Gile
(2008).  As such, it is referred to as the Sequential Sampling - Population Size Estimate (SS-PSE).
It uses the order of selection of the sample to provide information
on the distribution of network sizes over the population members.
}
\details{
It computes the posterior predictive distribution for each reported network size 
and computes the percentile rank of the reported network size within that posterior.
The percentile rank should be about 0.5 for a well specified model, but could be close to
uniform if there is little information about the reported network size.
The percentile ranks should not be extreme (e.g., close to zero or one) on a consistent basis
as this indicates a misspecified model.
}
\examples{

\dontrun{
data(fauxmadrona)
# Here interval=1 so that it will run faster. It should be higher in a 
# real application.
fit <- posteriorsize(fauxmadrona, median.prior.size=1000,
                                 burnin=20, interval=1, samplesize=100)
summary(fit)
# Let's look at some MCMC diagnostics
pospreddeg(fit)
}

}
\references{
Gile, Krista J. (2008) \emph{Inference from Partially-Observed Network
Data}, Ph.D. Thesis, Department of Statistics, University of Washington.

Gile, Krista J. and Handcock, Mark S. (2010) \emph{Respondent-Driven
Sampling: An Assessment of Current Methodology}, Sociological Methodology
40, 285-327.

Gile, Krista J. and Handcock, Mark S. (2014) \pkg{sspse}: Estimating Hidden 
Population Size using Respondent Driven Sampling Data
R package, Los Angeles, CA.  Version 0.5, \url{https://hpmrg.org/sspse/}.

Handcock MS (2003).  \pkg{degreenet}: Models for Skewed Count Distributions
Relevant to Networks.  Statnet Project, Seattle, WA.  Version 1.2,
\url{https://statnet.org/}.

Handcock, Mark S., Gile, Krista J. and Mar, Corinne M. (2014)
\emph{Estimating Hidden Population Size using Respondent-Driven Sampling
Data}, Electronic Journal of Statistics, 8, 1, 1491-1521

Handcock, Mark S., Gile, Krista J. and Mar, Corinne M. (2015)
\emph{Estimating the Size of Populations at High Risk for HIV using Respondent-Driven 
Sampling Data}, Biometrics.
}
\seealso{
The model fitting function \code{\link{posteriorsize}},
\code{\link[graphics]{plot}}.
}
