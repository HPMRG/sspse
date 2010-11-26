\name{size-package}
\alias{size-package}
\docType{package}
\title{
Estimating Hidden Population Size using Respondent Driven Sampling Data
}
\description{
\pkg{size} is a collection of functions to 
estimate population size and disease prevalence using 
respondent-driven sampling data. 
For a list of functions type: help(package='size')

For a complete list of the functions, use \code{library(help="size")} 
or read the rest of the manual.
% For a simple demonstration, use \code{demo(packages="size")}.

When publishing results obtained using this package the original authors
are to be cited as:

Gile, Krista G. 2003 
\emph{size: Software tools Estimating Hidden Population Size}
\cr
\url{http://statnet.org}.

All programs derived from this package must cite it.
For complete citation information, use\cr
\code{citation(package="size")}.
}
\details{
%Recent advances in the statistical modeling of random networks have had an
%impact on the empirical study of social networks. Statistical exponential
%family models (Strauss and Ikeda 1990) are a generalization of the Markov
%random network models introduced by Frank and Strauss (1986), which in turn
%derived from developments in spatial statistics (Besag, 1974). These models
%recognize the complex dependencies within relational data structures.  To
%date, the use of stochastic network models for networks has been limited by
%three interrelated factors: the complexity of realistic models, the lack of
%simulation tools for inference and validation, and a poor understanding of the
%inferential properties of nontrivial models.
%
%This manual introduces software tools for the representation, visualization,
%and analysis of network data that address each of these previous shortcomings.
%The package relies on the \code{\link[network]{network}} package which allows networks to be
%represented in \R. The \pkg{size} package allows maximum likelihood estimates of
%exponential random network models to be calculated using Markov Chain Monte
%Carlo.  The package also provides tools for plotting networks, simulating
%networks and assessing model goodness-of-fit.

For detailed information on how to download and install the software,
go to the \pkg{size} website:
\url{http://statnet.org}.
A tutorial, support newsgroup, references and links to further resources are provided there.
}
\author{
Krista J. Gile \email{krista.gile@nuffield.ox.ac.uk},\cr
Mark S. Handcock \email{handcock@stat.washington.edu}

Maintainer: Krista J. Gile \email{krista.gile@nuffield.ox.ac.uk}
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

\keyword{ package }
\keyword{ models }
