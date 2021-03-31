#' Estimating hidden population size using RDS data
#' 
#' \code{\link{posteriorsize}} computes the posterior distribution of the
#' population size based on data collected by Respondent Driven Sampling. The
#' approach approximates the RDS via the Sequential Sampling model of Gile
#' (2008). As such, it is referred to as the Sequential Sampling - Population Size Estimate (SS-PSE).
#' It uses the order of selection of the sample to provide information
#' on the distribution of network sizes over the population members.
#' 
#' @param s either a vector of integers or an \code{rds.data.frame} providing network 
#' size information.
#' If a \code{rds.data.frame} is passed and \code{visibility=TRUE}, the default, then
#' the measurement error model is to used, whereby latent visibilities are used in place 
#' of the reported network sizes as the size variable. If a vector of integers is passed these 
#' are the network sizes in sequential order of recording (and the measurement model is not used).
#' @param s2 either a vector of integers or an \code{rds.data.frame} providing network 
#' size information for a second RDS sample subsequent to the first RDS recorded in \eqn{s}.
#' If a \code{rds.data.frame} is passed and \code{visibility=TRUE}, the default, then
#' the measurement error model is to used, whereby latent visibilities are used in place 
#' of the reported network sizes as the size variable. If a vector of integers is passed these 
#' are the network sizes in sequential order of recording (and the measurement model is not used).
#' @param previous character; optionally, the name of the variable in \eqn{s2}
#' indicating if the corresponding unit was sampled in the first RDS.
#' @param median.prior.size scalar; A hyperparameter being the mode of the
#' prior distribution on the population size.
#' @param interval count; the number of proposals between sampled statistics.
#' @param burnin count; the number of proposals before any MCMC sampling is
#' done. It typically is set to a fairly large number.
#' @param maxN integer; maximum possible population size. By default this is
#' determined from an upper quantile of the prior distribution.
#' @param K count; the maximum visibility for an individual. This is usually
#' calculated as \code{round(stats::quantile(s,0.80))}. It applies to network sizes and (latent) visibilities.
#' If logical and FALSE then the K is unbounded but set to compute the visibilities.
#' @param samplesize count; the number of Monte-Carlo samples to draw to
#' compute the posterior. This is the number returned by the
#' Metropolis-Hastings algorithm.The default is 1000.
#' @param quartiles.prior.size vector of length 2; A pair of hyperparameters
#' being the lower and upper quartiles of the prior distribution on the
#' population size. For example, \cr \code{quartiles.prior.size=c(1000,4000)}
#' corresponds to a prior where the lower quartile (25\%) is 1000 and the upper
#' (75\%) is 4000.
#' @param mean.prior.size scalar; A hyperparameter being the mean of the prior
#' distribution on the population size.
#' @param mode.prior.size scalar; A hyperparameter being the mode of the prior
#' distribution on the population size.
#' @param priorsizedistribution character; the type of parametric distribution
#' to use for the prior on population size. The options are \code{beta} (for a
#' Beta prior on the sample proportion (i.e. \eqn{n/N})), \code{flat}
#' (uniform), \code{nbinom} (Negative-Binomial), and \code{pln}
#' (Poisson-log-normal). The default is \code{beta}.
#' @param effective.prior.df scalar; A hyperparameter being the effective
#' number of samples worth of information represented in the prior distribution
#' on the population size. By default this is 1, but it can be greater (or
#' less!) to allow for different levels of uncertainty.
#' @param sd.prior.size scalar; A hyperparameter being the standard deviation
#' of the prior distribution on the population size.
#' @param mode.prior.sample.proportion scalar; A hyperparameter being the mode
#' of the prior distribution on the sample proportion \eqn{n/N}.
#' @param alpha scalar; A hyperparameter being the first parameter of the beta
#' prior model for the sample proportion. By default this is NULL, meaning that
#' 1 is chosen. it can be any value at least 1 to allow for different levels of
#' uncertainty.
#' @param visibilitydistribution count; the parametric distribution to use for the
#' individual network sizes (i.e., degrees). The options are \code{cmp},
#' \code{nbinom}, and \code{pln}.  These correspond to the
#' Conway-Maxwell-Poisson, Negative-Binomial, and Poisson-log-normal. The
#' default is \code{cmp}.
#' @param mean.prior.visibility scalar; A hyper parameter being the mean visibility for
#' the prior distribution for a randomly chosen person. The prior has this
#' mean.
#' @param sd.prior.visibility scalar; A hyper parameter being the standard deviation
#' of the visibility for a randomly chosen person.  The prior has this
#' standard deviation.
#' @param max.sd.prior.visibility scalar; The maximum allowed value of \code{sd.prior.visibility}.
#' If the passed or computed value is higher, it is reduced to this value.
#' This is done for numerical stability reasons.
#' @param df.mean.prior.visibility scalar; A hyper parameter being the degrees-of-freedom
#' of the prior for the mean. This gives the equivalent sample size that would
#' contain the same amount of information inherent in the prior.
#' @param df.sd.prior.visibility scalar; A hyper parameter being the degrees-of-freedom of
#' the prior for the standard deviation. This gives the equivalent sample size
#' that would contain the same amount of information inherent in the prior for
#' the standard deviation.
#' @param beta0.mean.prior scalar; A hyper parameter being the mean of the 
#' beta0 parameter distribution in the model for the number of recruits.
#' @param beta1.mean.prior scalar; A hyper parameter being the mean of the 
#' beta1 parameter distribution in the model for the number of recruits.
#' @param beta0.sd.prior scalar; A hyper parameter being the standard deviation of the 
#' beta0 parameter distribution in the model for the number of recruits.
#' @param beta1.sd.prior scalar; A hyper parameter being the standard deviation of the 
#' beta0 parameter distribution in the model for the number of recruits.
#' @param mem.optimism.prior scalar; A hyper parameter being the mean of the 
#' distribution of the optimism parameter.
#' @param df.mem.optimism.prior scalar; A hyper parameter being the degrees-of-freedom
#' of the prior for the optimism parameter. This gives the equivalent sample size that would
#' contain the same amount of information inherent in the prior.
#' @param mem.scale.prior scalar; A hyper parameter being the scale of the concentration of
#' baseline negative binomial measurement error model.
#' @param df.mem.scale.prior scalar; A hyper parameter being the degrees-of-freedom of
#' the prior for the standard deviation of the dispersion parameter in the visibility model.
#' This gives the equivalent sample size
#' that would contain the same amount of information inherent in the prior for
#' the standard deviation.
#' @param mem.overdispersion scalar; A parameter being the overdispersion of the negative binomial
#' distribution that is the baseline for the measurement error model.
#' @param visibility logical; Indicate if the measurement error model
#' is to be used, whereby latent visibilities are used in place of the reported 
#' network sizes as the unit size variable. If \code{TRUE} then a \code{rds.data.frame}
#' need to be passed to provide the RDS information needed for the measurement error model.
#' @param type.impute The type of imputation to use for the summary visibilities 
#' (returned in the component \code{visibilities}. The imputes are based on the posterior 
#' draws of the visibilities. 
#' It can be of type \code{distribution}, \code{mode},\code{median}, or \code{mean} 
#' with \code{median} the default, being the posterior median of the visibility for that person.
#' @param Np integer; The overall visibility distribution is a mixture of the
#' \code{Np} rates for \code{1:Np} and a parametric visibility distribution model
#' truncated below \code{Np}. Thus the model fits the proportions of the
#' population with visibility \code{1:Np} each with a separate parameter. This
#' should adjust for an lack-of-fit of the parametric visibility distribution model
#' at lower visibilities, although it also changes the model away from the
#' parametric visibility distribution model.
#' @param n integer; the number of people in the sample. This is usually computed from
#' \eqn{s} automatically and not usually specified by the user.
#' @param n2 integer; If \eqn{s2} is specified, this is the number of people in the second sample. 
#' This is usually computed from
#' \eqn{s} automatically and not usually specified by the user.
#' @param muproposal scalar; The standard deviation of the proposal
#' distribution for the mean visibility.
#' @param nuproposal scalar; The standard deviation of the proposal
#' distribution for the CMP scale parameter that determines the standard deviation of the visibility.
#' @param beta0proposal scalar; The standard deviation of the proposal
#' distribution for the beta0 parameter of the recruit model.
#' @param beta1proposal scalar; The standard deviation of the proposal
#' distribution for the beta1 parameter of the recruit model.
#' @param memmuproposal scalar; The standard deviation of the proposal
#' distribution for the log of the optimism parameter (that is, gamma).
#' @param memscaleproposal scalar; The standard deviation of the proposal
#' distribution for the log of the s.d. in the optimism model.
#' @param burnintheta count; the number of proposals in the Metropolis-Hastings
#' sub-step for the visibility distribution parameters (\eqn{\theta}) before any
#' MCMC sampling is done. It typically is set to a modestly large number.
#' @param burninbeta count; the number of proposals in the Metropolis-Hastings
#' sub-step for the visibility distribution parameters (\eqn{\beta}) before any
#' MCMC sampling is done. It typically is set to a modestly large number.
#' @param parallel count; the number of parallel processes to run for the
#' Monte-Carlo sample.  This uses MPI or PSOCK. The default is 1, that is not to
#' use parallel processing.
#' @param parallel.type The type of parallel processing to use. The options are
#' "PSOCK" or "MPI". This requires the corresponding type to be installed.
#' The default is "PSOCK".
#' @param seed integer; random number integer seed.  Defaults to \code{NULL} to
#' use whatever the state of the random number generator is at the time of the
#' call.
#' @param maxbeta scalar; The maximum allowed value of the \code{beta} parameter.
#' If the implied or computed value is higher, it is reduced to this value.
#' This is done for numerical stability reasons.
#' @param supplied list; If supplied, is a list with components \code{maxN} and
#' \code{sample}. In this case \code{supplied} is a matrix with a column named
#' \code{N} being a sample from a prior distribution for the population size.
#' The value \code{maxN} specifies the maximum value of the population size, a
#' priori.
#' @param max.coupons The number of recruitment coupons distributed to each 
#' enrolled subject (i.e. the maximum number of recruitees for any subject).
#' By default it is taken by the attribute or data, else the maximum recorded number of coupons.
#' @param recruit.time vector; An optional value for the data/time that the person was interviewed.
#' It needs to resolve as a numeric vector with number of elements the number
#' of rows of the data with non-missing values of the network variable. If it
#' is a character name of a variable in the data then that variable is used.
#' If it is NULL then the sequence number of the recruit in the data is used.
#' If it is NA then the recruitment is not used in the model.
#' Otherwise, the recruitment time is used in the model to better predict
#' the visibility of the person.
#' @param recruit.time2 vector; An optional value for the data/time that the person in the second RDS survey was interviewed.
#' It needs to resolve as a numeric vector with number of elements the number
#' of rows of the data with non-missing values of the network variable. If it
#' is a character name of a variable in the data then that variable is used.
#' If it is NULL, the default, then the sequence number of the recruit in the data is used.
#' If it is NA then the recruitment is not used in the model.
#' Otherwise, the recruitment time is used in the model to better predict
#' the visibility of the person.
#' @param include.tree logical; If \code{TRUE}, 
#' augment the reported network size by the number of recruits and one for the recruiter (if any).
#' This reflects a more accurate value for the visibility, but is not the self-reported degree.
#' In particular, it typically produces a positive visibility (compared to a
#' possibility zero self-reported degree). 
#' @param unit.scale numeric; If not \code{NULL} it sets the numeric value of the scale parameter
#' of the distribution of the unit sizes.
#' For the negative binomial, it is the multiplier on the variance of the negative binomial 
#' compared to a Poisson (via the Poisson-Gamma mixture representation). Sometimes the scale is 
#' unnaturally large (e.g. 40) so this give the option of fixing it (rather than using
#' the MLE of it). The model is fit with the parameter fixed at this passed value.
#' @param optimism logical; If \code{TRUE} then add a term to the model allowing
#' the (proportional) inflation of the self-reported degrees relative to the unit sizes.
#' @param reflect.time logical; If \code{TRUE} then the \code{recruit.time} is the time before the 
#' end of the study (instead of the time since the survey started or chronological time).
#' @param equalize logical; If \code{TRUE} and the capture-recapture model is used, adjusts for gross differences in the 
#' reported network sizes between the two samples.
#' @param verbose logical; if this is \code{TRUE}, the program will print out
#' additional information, including goodness of fit statistics.
#' @return \code{\link{posteriorsize}} returns a list consisting of the
#' following elements:
#'\item{pop}{vector; The final posterior draw for the
#' degrees of the population. The first \eqn{n} are the sample in sequence and
#' the reminder are non-sequenced.}
#'\item{K}{count; the maximum visibility for an
#' individual. This is usually calculated as twice the maximum observed
#' degree.}
#'\item{n}{count; the sample size.}
#'\item{samplesize}{count; the
#' number of Monte-Carlo samples to draw to compute the posterior. This is the
#' number returned by the Metropolis-Hastings algorithm.The default is 1000.}
#'\item{burnin}{count; the number of proposals before any MCMC sampling is
#' done. It typically is set to a fairly large number.}
#'\item{interval}{count; the number of proposals between sampled statistics.}
#'\item{mu}{scalar; The
#' hyper parameter \code{mean.prior.visibility} being the mean visibility for the prior
#' distribution for a randomly chosen person. The prior has this mean.}
#'\item{sigma}{scalar; The hyper parameter \code{sigma} being the
#' standard deviation of the visibility for a randomly chosen person. The prior has
#' this standard deviation.}
#'\item{df.mean.prior.visibility}{scalar; A hyper parameter
#' being the degrees-of-freedom of the prior for the mean. This gives the
#' equivalent sample size that would contain the same amount of information
#' inherent in the prior.}
#'\item{df.sd.prior.visibility}{scalar; A hyper parameter being
#' the degrees-of-freedom of the prior for the standard deviation. This gives
#' the equivalent sample size that would contain the same amount of information
#' inherent in the prior for the standard deviation.}
#'\item{Np}{integer; The
#' overall visibility distribution is a mixture of the \code{1:Np} rates and a
#' parametric visibility distribution model truncated below Np. Thus the model fits
#' the proportions of the population with visibility \code{1:Np} each with a
#' separate parameter. This should adjust for an lack-of-fit of the parametric
#' visibility distribution model at lower visibilities, although it also changes the
#' model away from the parametric visibility distribution model.}
#'\item{muproposal}{scalar; The standard deviation of the proposal
#' distribution for the mean visibility.}
#'\item{nuproposal}{scalar; The standard
#' deviation of the proposal distribution for the CMP scale parameter of the
#' visibility distribution.}
#'\item{N}{vector of length 5; summary statistics for the posterior
#' population size.
#' \describe{
#'\item{MAP}{maximum aposteriori value of N}
#'\item{Mean AP}{mean aposteriori value of N}
#'\item{Median AP}{median aposteriori value of N}
#'\item{P025}{the 2.5th
#' percentile of the (posterior) distribution for the N. That is, the lower
#' point on a 95\% probability interval.}
#'\item{P975}{the 97.5th
#' percentile of the (posterior) distribution for the N. That is, the upper
#' point on a 95\% probability interval.} } }
#'\item{maxN}{integer; maximum
#' possible population size. By default this is determined from an upper
#' quantile of the prior distribution.}
#'\item{sample}{matrix of dimension
#' \code{samplesize}\eqn{\times} \code{10} matrix of summary statistics from
#' the posterior. This is also an object of class \code{mcmc} so it can be
#' plotted and summarized via the \code{mcmc.diagnostics} function in the
#' \code{ergm} package (and also the \code{coda} package). The statistics are:
#'\describe{
#'\item{N}{population size.}
#'\item{mu}{scalar; The mean
#' visibility for the prior distribution for a randomly chosen person. The prior
#' has this mean.}
#'\item{sigma}{scalar; The standard deviation of the visibility
#' for a randomly chosen person. The prior has this standard deviation.}
#'\item{visibility1}{scalar; the number of nodes of visibility 1 in the population (it
#' is assumed all nodes have visibility 1 or more).}
#'\item{lambda}{scalar; This is
#' only present for the \code{cmp} model. It is the \eqn{\lambda} parameter in
#' the standard parameterization of the Conway-Maxwell-Poisson model for the
#' visibility distribution.}
#'\item{nu}{scalar; This is only present for the
#' \code{cmp} model. It is the \eqn{\nu} parameter in the standard
#' parameterization of the Conway-Maxwell-Poisson model for the visibility
#' distribution.} } }
#'\item{sample}{matrix of dimension \code{samplesize}\eqn{\times} \code{n} matrix of 
#' posterior.draws from the unit size distribution for those in the survey.
#' The sample for the \code{i}th person is the \code{i}th column.}
#'\item{lpriorm}{vector; the vector of (log) prior
#' probabilities on each value of \eqn{m=N-n} - that is, the number of
#' unobserved members of the population. The values are
#' \code{n:(length(lpriorm)-1+n)}.}
#'\item{burnintheta}{count; the number of
#' proposals in the Metropolis-Hastings sub-step for the visibility distribution
#' parameters (\eqn{\theta}) before any MCMC sampling is done. It typically is
#' set to a modestly large number.}
#'\item{verbose}{logical; if this is
#' \code{TRUE}, the program printed out additional information, including
#' goodness of fit statistics.}
#'\item{predictive.visibility.count}{vector; a vector
#' of length the maximum visibility (\code{K}) (by default \cr \code{K=2*max(sample
#' visibility)}).  The \code{k}th entry is the posterior predictive number persons
#' with visibility \code{k}.  That is, it is the posterior predictive distribution
#' of the number of people with each visibility in the population.}
#'\item{predictive.visibility}{vector; a vector of length the maximum visibility
#' (\code{K}) (by default \cr \code{K=2*max(sample visibility)}).  The \code{k}th entry
#' is the posterior predictive proportion of persons with visibility \code{k}.
#' That is, it is the posterior predictive distribution of the proportion of
#' people with each visibility in the population.}
#'\item{MAP}{vector of length 6
#' of MAP estimates corresponding to the output \code{sample}. These are:
#'\describe{
#'\item{N}{population size.}
#'\item{mu}{scalar; The mean
#' visibility for the prior distribution for a randomly chosen person. The prior
#' has this mean.}
#'\item{sigma}{scalar; The standard deviation of the visibility
#' for a randomly chosen person. The prior has this standard deviation.}
#'\item{visibility1}{scalar; the number of nodes of visibility 1 in the population (it
#' is assumed all nodes have visibility 1 or more).}
#'\item{lambda}{scalar; This is
#' only present for the \code{cmp} model. It is the \eqn{\lambda} parameter in
#' the standard parameterization of the Conway-Maxwell-Poisson model for the
#' visibility distribution.}
#'\item{nu}{scalar; This is only present for the
#' \code{cmp} model. It is the \eqn{\nu} parameter in the standard
#' parameterization of the Conway-Maxwell-Poisson model for the visibility
#' distribution.} } }
#'\item{mode.prior.sample.proportion}{scalar; A
#' hyperparameter being the mode of the prior distribution on the sample
#' proportion \eqn{n/N}.}
#'\item{median.prior.size}{scalar; A hyperparameter
#' being the mode of the prior distribution on the population size.}
#'\item{mode.prior.size}{scalar; A hyperparameter being the mode of the prior
#' distribution on the population size.}
#'\item{mean.prior.size}{scalar; A
#' hyperparameter being the mean of the prior distribution on the population
#' size.}
#'\item{quartiles.prior.size}{vector of length 2; A pair of
#' hyperparameters being the lower and upper quartiles of the prior
#' distribution on the population size.}
#'\item{visibilitydistribution}{count; the
#' parametric distribution to use for the individual network sizes (i.e.,
#' visibilities). The options are \code{cmp}, \code{nbinom}, and \code{pln}.  These
#' correspond to the Conway-Maxwell-Poisson, Negative-Binomial, and
#' Poisson-log-normal. The default is \code{cmp}.}
#' \item{priorsizedistribution}{character; the type of parametric distribution
#' to use for the prior on population size. The options are \code{beta} (for a
#' Beta prior on the sample proportion (i.e. \eqn{n/N}), \code{nbinom}
#' (Negative-Binomial), \code{pln} (Poisson-log-normal), \code{flat} (uniform),
#' and \code{continuous} (the continuous version of the Beta prior on the
#' sample proportion. The default is \code{beta}. }
#' @section Details on priors: The best way to specify the prior is via the
#' hyperparameter \code{mode.prior.size} which specifies the mode of the prior
#' distribution on the population size. You can alternatively specify the
#' hyperparameter \code{median.prior.size} which specifies the median of the
#' prior distribution on the population size, or \code{mean.prior.sample
#' proportion} which specifies the mean of the prior distribution on the
#' proportion of the population size in the sample or \code{mode.prior.sample
#' proportion} which specifies the mode of the prior distribution on the
#' proportion of the population size in the sample. Finally, you can specify
#' \code{quartiles.prior.size} as a vector of length 2 being the pair of lower
#' and upper quartiles of the prior distribution on the population size.
#' @seealso network, statnet, degreenet
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
#' R package, Los Angeles, CA.  Version 0.5, \url{http://hpmrg.org}.
#' 
#' Handcock MS (2003).  \pkg{degreenet}: Models for Skewed Count Distributions
#' Relevant to Networks.  Statnet Project, Seattle, WA.  Version 1.2,
#' \url{http://statnetproject.org}.
#' 
#' Handcock, Mark S., Gile, Krista J. and Mar, Corinne M. (2014)
#' \emph{Estimating Hidden Population Size using Respondent-Driven Sampling
#' Data}, Electronic Journal of Statistics, 8, 1, 1491-1521
#' 
#' Handcock, Mark S., Gile, Krista J. and Mar, Corinne M. (2015)
#' \emph{Estimating the Size of Populations at High Risk for HIV using Respondent-Driven 
#' Sampling Data}, Biometrics.
#' @keywords models
#' @examples
#' data(fauxmadrona)
#' # Here interval=1 so that it will run faster. It should be higher in a 
#' # real application.
#' fit <- posteriorsize(fauxmadrona, median.prior.size=1000,
#'                                  burnin=100, interval=1, samplesize=100)
#' summary(fit)
#' @export posteriorsize
posteriorsize<-function(s,
                  s2=NULL, previous=NULL,
                  median.prior.size=NULL,
                  interval=10,
                  burnin=5000,
                  maxN=NULL,
                  K=FALSE,
                  samplesize=1000,
                  quartiles.prior.size=NULL,
                  mean.prior.size=NULL,
                  mode.prior.size=NULL,
                  priorsizedistribution=c("beta","flat","nbinom","pln","supplied"),
                  effective.prior.df=1,
                  sd.prior.size=NULL,
                  mode.prior.sample.proportion=NULL,
                  alpha=NULL,
                  visibilitydistribution=c("cmp","nbinom","pln"),
                  mean.prior.visibility=NULL, sd.prior.visibility=NULL, max.sd.prior.visibility=4,
                  df.mean.prior.visibility=1,df.sd.prior.visibility=3,
                  beta0.mean.prior=-3, beta1.mean.prior=0,
                  beta0.sd.prior=10, beta1.sd.prior=10,
                  mem.optimism.prior=NULL, df.mem.optimism.prior=5, 
                  mem.scale.prior=2, df.mem.scale.prior=10,
		  mem.overdispersion=15,
                  visibility=TRUE,
                  type.impute = c("median","distribution","mode","mean"),
                  Np=0,
                  n=NULL,
                  n2=NULL,
                  muproposal=0.1, 
                  nuproposal=0.15, 
                  beta0proposal=0.2, beta1proposal=0.001,
                  memmuproposal=0.1, memscaleproposal=0.15,
                  burnintheta=500,
                  burninbeta=50,
                  parallel=1, parallel.type="PSOCK", seed=NULL, 
                  maxbeta=90, 
                  supplied=list(maxN=maxN),
                  max.coupons=NULL,
                  recruit.time=NULL,recruit.time2=NULL,
                  include.tree=TRUE, unit.scale=FALSE, 
                  optimism = TRUE,
                  reflect.time=FALSE,
                  equalize=TRUE,
                  verbose=TRUE){
#
  visibilitydistribution=match.arg(visibilitydistribution)
  posfn <- switch(visibilitydistribution,
                  nbinom=posnbinom,
                  pln=pospln,
                  cmp=poscmpwp,
                  poscmpwp)
  K.fixed <- K
  # If the passed "s" is an rds.rata.frame, extract out the components
  rc <- NULL
  nr <- 1
  recruit.times <- 1
  remvalues <- TRUE
  if(!methods::is(s,"rds.data.frame")){
   # a sequence is passed
   rds.data <- NULL
   visibility <- FALSE
   if(is.null(K.fixed)) K=max(s,na.rm=TRUE)
   if(is.null(n)) n=length(s)
  }else{
  # an rds.data.frame is passed
  rds.data <- s
  n <- nrow(rds.data)
  if(is.null(attr(rds.data,"network.size.variable")))
    stop("rds.data must have a network.size attribute.")
  nr <- RDS::get.number.of.recruits(rds.data)
  nw <- RDS::get.wave(rds.data)
  ns <- RDS::get.seed.id(rds.data)
  is.seed <- (RDS::get.rid(rds.data)=="seed")
  
  if(is.null(max.coupons)){
    max.coupons <- attr(rds.data,"max.coupons")
    if(is.null(max.coupons)){
      max.coupons <- max(nr,na.rm=TRUE)
    }
  }
  if(length(recruit.time)==1){
    if(is.character(recruit.time)){
      if(recruit.time=="wave"){
        recruit.times <- nw
      }else{
       recruit.times <- rds.data[[recruit.time]]
       if(methods::is(recruit.times,"POSIXt") | methods::is(recruit.times,"Date")){
        recruit.times <- as.numeric(recruit.times) / (24*60*60)
       }else{
        recruit.times <- as.numeric(recruit.times)
       }
      }
      recruit.time <- TRUE
    }else{
      if(is.na(recruit.time)){
        recruit.times <- rep(0,n)
        recruit.time <- FALSE
      }else{
        stop("The recruitment time should be a variable in the RDS data, or 'wave' to indicate the wave number or NA/NULL to indicate that the recruitment time is not available and/or used.")
      }
    }
  }else{
    if(length(recruit.time)==0 & is.null(recruit.time)){
      recruit.time <- 1:n
    }else{
      if(length(recruit.time)!=n | (!is.numeric(recruit.time) & !methods::is(recruit.time,"POSIXt") & !methods::is(recruit.time,"Date"))){
        stop("The recruitment time should be a variable in the RDS data, or 'wave' to indicate the wave number or NA/NULL to indicate that the recruitment time is not available and/or used.")
      }
    }
    if(length(recruit.time)==n & (methods::is(recruit.time,"POSIXt") | methods::is(recruit.time,"Date"))){
      recruit.times <- as.numeric(recruit.time) / (24*60*60)
    }else{
      recruit.times <- recruit.time
    }
    recruit.time <- TRUE
  }
  if(any(is.na(recruit.times))){
    med.index <- cbind(c(2,1:(n-1)),c(3,3:n,n))
    moving.median=function(i){stats::median(recruit.times[med.index[i,]],na.rm=TRUE)}
    while(any(is.na(recruit.times))){
      for(i in which(is.na(recruit.times))){recruit.times[i] <- moving.median(i)}
    }
  }
# gap <- diff(sort(recruit.times))
# gap <- min(gap[gap > 0])
# recruit.times <- recruit.times + 0.01*(1:n)*gap/(n+1)
  recruit.times <- recruit.times - min(recruit.times)
  if(reflect.time){
    recruit.times <- max(recruit.times)-recruit.times
  }
  network.size <- as.numeric(rds.data[[attr(rds.data,"network.size.variable")]])
  remvalues <- is.na(network.size)
  if(any(remvalues)){
    warning(paste(sum(remvalues),"of",nrow(rds.data),
                  "network sizes were missing. These will be imputed from the marginal distribution"), call. = FALSE)
  }
  
  if(!is.null(K) & is.logical(K) & (K==FALSE)){
   if(visibility){
#   K.fixed <- max(network.size[!remvalues])
#   K.fixed <- round(stats::quantile(network.size[!remvalues],0.99))
    K.fixed <- round(stats::quantile(network.size[!remvalues],0.95))
   }else{
    K.fixed <- NULL
   }
  }
  if(is.null(K.fixed)){
    if(length(network.size[!remvalues])>0){
      K <- round(stats::quantile(network.size[!remvalues],0.95))
#     K <- round(stats::quantile(network.size[!remvalues],0.99))
    }
  }
  
  #Augment the reported network size by the number of recruits and the recruiter (if any).
  if(include.tree){
    nsize <- pmax(network.size,nr+!is.seed)
  }else{
    nsize <- network.size
  }
  
  gmean <- HT.estimate(RDS::vh.weights(nsize[!is.na(nsize)]),nsize[!is.na(nsize)])
  if(is.na(gmean)) gmean <- 38
  
  recruit.times.order <- order(order(recruit.times))
  recruit.times.order.notrem <- order(order(recruit.times[!remvalues]))
  s <- nsize[order(recruit.times)]
  recruit.times <- recruit.times[order(recruit.times)]
  nr <- nr[order(recruit.times)]
  }
  # End of measurement model information extraction for the first RDS

  # If the passed "s2" is an rds.rata.frame, extract out the components
  rc2 <- NULL
  nr2 <- 1
  recruit.times2 <- 1
  remvalues2 <- TRUE
  rds.data2 <- NULL
  if(!is.null(s2)){
  if(!methods::is(s2,"rds.data.frame")){
   # a sequence is passed
   visibility <- FALSE
   if(is.null(K.fixed)) K=max(c(s,s2),na.rm=TRUE)
   if(is.null(n2)) n2=length(s2)
  }else{
   # an rds.data.frame is passed
   rds.data2 <- s2
   n2 <- nrow(rds.data2)
   if(is.null(attr(rds.data2,"network.size.variable"))){
     stop("The second rds.data must have a network.size attribute.")
   }
   nr2 <- RDS::get.number.of.recruits(rds.data2)
   nw2 <- RDS::get.wave(rds.data2)
   ns2 <- RDS::get.seed.id(rds.data2)
   is.seed2 <- (RDS::get.rid(rds.data2)=="seed")
   
   max.coupons2 <- attr(rds.data2,"max.coupons")
   if(is.null(max.coupons2)){
     max.coupons2 <- max(nr2,na.rm=TRUE)
   }
   if(length(recruit.time2)==1){
     if(is.character(recruit.time2)){
       if(recruit.time2=="wave"){
         recruit.times2 <- nw2
       }else{
        recruit.times2 <- rds.data2[[recruit.time2]]
        if(methods::is(recruit.times2,"POSIXt") | methods::is(recruit.times2,"Date")){
         recruit.times2 <- as.numeric(recruit.times2) / (24*60*60)
        }else{
         recruit.times2 <- as.numeric(recruit.times2)
        }
       }
       recruit.time2 <- TRUE
     }else{
       if(is.na(recruit.time2)){
         recruit.times2 <- rep(0,n2)
         recruit.time2 <- FALSE
       }else{
         stop("The recruitment time should be a variable in the second RDS data, or 'wave' to indicate the wave number or NA/NULL to indicate that the recruitment time is not available and/or used.")
       }
     }
   }else{
     if(length(recruit.time2)==0 & is.null(recruit.time2)){
       recruit.time2 <- 1:n2
     }else{
       if(length(recruit.time2)!=n2 | (!is.numeric(recruit.time2) & !methods::is(recruit.time2,"POSIXt") & !methods::is(recruit.time2,"Date"))){
         stop("The recruitment time should be a variable in the second RDS data, or 'wave' to indicate the wave number or NA/NULL to indicate that the recruitment time is not available and/or used.")
       }
     }
     if(length(recruit.time2)==n & (methods::is(recruit.time2,"POSIXt") | methods::is(recruit.time2,"Date"))){
       recruit.times2 <- as.numeric(recruit.time2) / (24*60*60)
     }else{
       recruit.times2 <- recruit.time2
     }
     recruit.time2 <- TRUE
   }
   if(any(is.na(recruit.times2))){
     med.index <- cbind(c(2,1:(n2-1)),c(3,3:n2,n2))
     moving.median=function(i){stats::median(recruit.times2[med.index[i,]],na.rm=TRUE)}
     while(any(is.na(recruit.times2))){
       for(i in which(is.na(recruit.times2))){recruit.times2[i] <- moving.median(i)}
     }
   }
 # gap <- diff(sort(recruit.times))
 # gap <- min(gap[gap > 0])
 # recruit.times <- recruit.times + 0.01*(1:n)*gap/(n+1)
   recruit.times2 <- recruit.times2 - min(recruit.times2)
   if(reflect.time){
     recruit.times2 <- max(recruit.times2)-recruit.times2
   }
   network.size2 <- as.numeric(rds.data2[[attr(rds.data2,"network.size.variable")]])
   remvalues2 <- is.na(network.size2)
   if(any(remvalues2)){
     warning(paste(sum(remvalues2),"of",nrow(rds.data2),
                   "network sizes were missing in the second RDS data set. These will be imputed from the marginal distribution"), call. = FALSE)
   }
   
   if(equalize){
    if(sum(!remvalues) >= sum(!remvalues2)){
     a <- rank(network.size2[!remvalues2],ties.method="random")
     network.size2[!remvalues2] <- sort((network.size[!remvalues])[1:sum(!remvalues2)])[a]
    }else{
     a <- round(sum(!remvalues)*rank(network.size2[!remvalues2],ties.method="random")/sum(!remvalues2))
     network.size2[!remvalues2] <- sort(network.size[!remvalues])[a]
    }
    message(sprintf("Adjusting for the gross differences in the reported network sizes between the two samples.\n"),appendLF=FALSE)
   }

   if(!is.null(K) & is.logical(K) & (K==FALSE)){
    if(visibility){
#    K.fixed <- max(network.size2[!remvalues2])
     K.fixed <- max(K.fixed, round(stats::quantile(network.size2[!remvalues2],0.99)))
    }else{
     K.fixed <- NULL
    }
   }
   if(is.null(K.fixed) & length(network.size2[!remvalues2])>0){
     K <- max(K, round(stats::quantile(network.size2[!remvalues2],0.95)))
   }
   
   #Augment the reported network size by the number of recruits and the recruiter (if any).
   if(include.tree){
     nsize2 <- pmax(network.size2,nr2+!is.seed2)
   }else{
     nsize2 <- network.size2
   }
   
   gmean <- HT.estimate(RDS::vh.weights(nsize2[!is.na(nsize2)]),nsize2[!is.na(nsize2)])
   if(is.na(gmean)) gmean <- 38
   
   recruit.times2.order <- order(order(recruit.times2))
   recruit.times2.order.notrem <- order(order(recruit.times2[!remvalues2]))
   s2 <- nsize2[order(recruit.times2)]
   recruit.times2 <- recruit.times2[order(recruit.times2)]
   nr2 <- nr2[order(recruit.times2)]
  }}
  # End of measurement model information extraction for the second survey
  # End of measurement model information extraction

  remvalues <- is.na(s)
  if(sum(!remvalues) < length(s)){
   warning(paste(length(s)-sum(!remvalues),"of",length(s),
          "sizes values were missing and were removed."), call. = FALSE)
   s <- s[!remvalues]
   n <- length(s)
  }
  s.prior <- s
  if(!is.null(s2)){
    if(is.null(rds.data2[[previous]])){
      stop("The argument 'previous' must have a variable in the second RDS data set indicating if the corresponding unit was sampled in the first RDS.")
    }
    rc <- rds.data2[[previous]]
    rc <- rc[order(recruit.times2)]
    if(!is.logical(rc) | length(s2)!=length(rc)){
      stop("The argument 'previous' must have a variable in the second RDS data set that is a logical vector of the same length as s2, indicating if the corresponding unit was sampled in the first RDS.")
    }
#   rc <- !(rc==0)
    remvalues2 <- is.na(s2)
    if(sum(!remvalues2) < length(s2)){
     warning(paste(length(s2)-sum(!remvalues2),"of",length(s2),
            "sizes values from the second RDS were missing and were removed."), call. = FALSE)
     s2 <- s2[!remvalues2]
     rc <- rc[!remvalues2]
     n2 <- length(s2)
    }
    s.prior <- c(s.prior,s2[!rc])
  }
  priorsizedistribution=match.arg(priorsizedistribution)
  # Extract the population size from the network if it is not set.
  if(is.null(mode.prior.sample.proportion)
   & is.null(median.prior.size)
   & is.null(mean.prior.size)
   & is.null(mode.prior.size)
   & methods::is(rds.data,"rds.data.frame")
     ){
   median.prior.size <- attr(rds.data, "population.size.mid")
   if(methods::is(rds.data2,"rds.data.frame")){
     if(is.null(median.prior.size)){
      median.prior.size <- attr(rds.data2, "population.size.mid")
     }else{
      median.prior.size2 <- attr(rds.data2, "population.size.mid")
      if(!is.null(median.prior.size2)){
       median.prior.size <- 0.5*(median.prior.size+median.prior.size2)
      }
     }
   }
    if(!is.null(median.prior.size)){
     warning(paste("The median of the prior distribution of the population size is set to", 
                    median.prior.size), call. = FALSE)
    }
  }
  if(priorsizedistribution=="nbinom" && missing(mean.prior.size)){
    stop("You need to specify 'mean.prior.size', and possibly 'sd.prior.size' if you use the 'nbinom' prior.") 
  }
  if(is.null(K.fixed) & visibility){
    message(sprintf("Initial cap on influence of the visibility is K = %d.\n",K),appendLF=FALSE)
    K=round(stats::quantile(s.prior,0.90))+1
    degs <- s.prior
    degs[degs>K] <- K
    degs[degs==0]<-1
    ds<-degs
    isnas <- is.na(degs)
    degs <- sum(!isnas)*(degs)/sum(degs,na.rm=TRUE)
    weights <- (1/degs)
    weights[is.na(weights)] <- 0
    mean.pd <- sum(ds*weights)/sum(weights)
    sd.pd <- sum(ds*ds*weights)/sum(weights)
    sd.pd <- sqrt(sd.pd - mean.pd^2)
    if(sd.pd > max.sd.prior.visibility*mean.pd){
     sd.pd <- min(max.sd.prior.visibility*mean.pd, sd.pd)
    }
    xv <- ds
#   xp <- weights*ds
    xp <- weights
    xp <- length(xp)*xp/sum(xp)
    txp <- tapply(xp,xv,sum)
    txv <- tapply(xv,xv,stats::median)
    fit <- cmpmle(txv,txp,cutoff=1,cutabove=K-1,guess=c(mean.pd, sd.pd))
    y=dcmp.natural(v=fit,x=(0:max(s.prior)))
    K=(0:max(s.prior))[which.max(cumsum(y)>0.99)]
#   K=round(stats::quantile(s,0.99))
  }
  if(!is.null(K.fixed) & abs(df.mem.optimism.prior-10001)>0.001){K <- K.fixed}
  if(visibility & is.null(mem.optimism.prior)){
    degs <- network.size
#   degs[degs>K.fixed] <- K.fixed
    degs[degs==0]<-1
    ds<-degs
    isnas <- is.na(degs)
    degs <- sum(!isnas)*(degs)/sum(degs,na.rm=TRUE)
#   weights <- (1/degs)
    weights <- degs-degs+1
    weights[is.na(weights)] <- 0
    mean.pd <- sum(ds*weights)/sum(weights)
#   sd.pd <- sum(ds*ds*weights)/sum(weights)
#   sd.pd <- sqrt(sd.pd - mean.pd^2)
#   if(sd.pd > max.sd.prior.visibility*mean.pd){
#    sd.pd <- min(max.sd.prior.visibility*mean.pd, sd.pd)
#   }
    message(sprintf("The mean of the recorded personal network sizes is %f.\n",mean.pd),appendLF=FALSE)
    xv <- ds
    xp <- weights
    xp <- length(xp)*xp/sum(xp)
    txp <- tapply(xp,xv,sum)
    txv <- tapply(xv,xv,stats::median)
    mem.optimism.prior <- nb0mle(xv=txv,xp=txp,cutabove=K.fixed-1)[1] / 8
#   mem.optimism.prior <- mean.pd / 8
    K <- 35
    df.mem.optimism.prior <- 10001
    message(sprintf("mem.optimism.prior set is to %f.\n",mem.optimism.prior),appendLF=FALSE)
  }
  if(visibility){
    message(sprintf("The cap on the visibility distribution is K = %d.\n",K),appendLF=FALSE)
  }else{
    message(sprintf("The cap on influence of the personal network size is K = %d.\n",K),appendLF=FALSE)
  }
  if(is.null(mean.prior.visibility)){
    degs <- s.prior
    degs[degs>K] <- K
    degs[degs==0]<-1
    ds<-degs
    isnas <- is.na(degs)
    degs <- sum(!isnas)*(degs)/sum(degs,na.rm=TRUE)
    weights <- (1/degs)
    weights[is.na(weights)] <- 0
    mean.prior.visibility <- sum(ds*weights)/sum(weights)
    if(is.null(sd.prior.visibility)){
     sd.prior.visibility <- sum(ds*ds*weights)/sum(weights)
     sd.prior.visibility <- sqrt(sd.prior.visibility - mean.prior.visibility^2)
    }
    xv <- ds
#   xp <- weights*ds
    xp <- weights
    xp[is.na(xp)] <- 0
    xp <- length(xp)*xp/sum(xp)
    txp <- tapply(xp,xv,sum)
    txv <- tapply(xv,xv,stats::median)
    fit <- cmpmle(txv,txp,cutoff=1,cutabove=K-1,
            guess=c(mean.prior.visibility,sd.prior.visibility))
    fit <- cmp.to.mu.sd(fit,max.mu=5*mean.prior.visibility,force=TRUE)
    if(verbose){
      message(sprintf("The preliminary empirical value of the mean of the prior distribution for visibility is %f.\n",mean.prior.visibility),appendLF=FALSE)
      message(sprintf("The preliminary empirical value of the s.d. of the prior distribution for visibility is %f.\n",sd.prior.visibility),appendLF=FALSE)
    }
    mean.prior.visibility = fit[1]
    sd.prior.visibility = fit[2]
  }else{
    if(is.null(sd.prior.visibility)){sd.prior.visibility <- sqrt(mean.prior.visibility)}
  }
  if(verbose){
    message(sprintf("The mean of the prior distribution for visibility is %f.\n",mean.prior.visibility),appendLF=FALSE)
    message(sprintf("The s.d. of the prior distribution for visibility is %f.\n",sd.prior.visibility),appendLF=FALSE)
  }
  if(sd.prior.visibility > max.sd.prior.visibility*mean.prior.visibility){
    sd.prior.visibility <- min(max.sd.prior.visibility*mean.prior.visibility, sd.prior.visibility)
    message(sprintf("The suggested s.d. of the prior distribution for visibility is too large and has been reduced to the more reasonable %f.\n",sd.prior.visibility),appendLF=FALSE)
  }
  ### are we running the job in parallel (parallel > 1), if not just 
  #   call the visibility specific function
  if(parallel==1){
      Cret <- posfn(s=s,s2=s2,rc=rc,K=K,maxN=maxN,
                    mean.prior.visibility=mean.prior.visibility,df.mean.prior.visibility=df.mean.prior.visibility,
                    sd.prior.visibility=sd.prior.visibility,df.sd.prior.visibility=df.sd.prior.visibility,
                    beta0.mean.prior=beta0.mean.prior, beta1.mean.prior=beta1.mean.prior,
                    beta0.sd.prior=beta0.sd.prior, beta1.sd.prior=beta1.sd.prior,
                    mem.optimism.prior=mem.optimism.prior, df.mem.optimism.prior=df.mem.optimism.prior,
                    mem.scale.prior=mem.scale.prior, df.mem.scale.prior=df.mem.scale.prior,
		    mem.overdispersion=mem.overdispersion,
                    muproposal=muproposal, nuproposal=nuproposal, 
                    beta0proposal=beta0proposal, beta1proposal=beta1proposal,
                    memmuproposal=memmuproposal, memscaleproposal=memscaleproposal,
                    visibility=visibility,
                    Np=Np,
                    samplesize=samplesize,burnin=burnin,interval=interval,
                    burnintheta=burnintheta,
                    burninbeta=burninbeta,
                    priorsizedistribution=priorsizedistribution,
                    mean.prior.size=mean.prior.size, sd.prior.size=sd.prior.size,
                    mode.prior.sample.proportion=mode.prior.sample.proportion,
                    median.prior.size=median.prior.size,
                    mode.prior.size=mode.prior.size,
                    quartiles.prior.size=quartiles.prior.size,
                    effective.prior.df=effective.prior.df,
                    alpha=alpha,
                    seed=seed,
                    supplied=supplied,
                    num.recruits=nr[!remvalues],
                    recruit.times=recruit.times[!remvalues],
                    num.recruits2=nr2[!remvalues2],
                    recruit.times2=recruit.times2[!remvalues2],
                    max.coupons=max.coupons,
                    maxbeta=maxbeta)
  }else{
  ### since running job in parallel, start vm (if not already running)
    cl <- beginparallel(parallel,type=parallel.type)
    ### divide the samplesize by the number of parallel runs (number of MCMC samples)
    samplesize.parallel=round(samplesize/parallel)
    ### cluster call, send following to each of the virtual machines, posnbinom function
    ### with it's arguments
    outlist <- parallel::clusterCall(cl, posfn,
      s=s,s2=s2,rc=rc,K=K,maxN=maxN,
      mean.prior.visibility=mean.prior.visibility,df.mean.prior.visibility=df.mean.prior.visibility,
      sd.prior.visibility=sd.prior.visibility,df.sd.prior.visibility=df.sd.prior.visibility,
      beta0.mean.prior=beta0.mean.prior, beta1.mean.prior=beta1.mean.prior,
      beta0.sd.prior=beta0.sd.prior, beta1.sd.prior=beta1.sd.prior,
      mem.optimism.prior=mem.optimism.prior, df.mem.optimism.prior=df.mem.optimism.prior,
      mem.scale.prior=mem.scale.prior, df.mem.scale.prior=df.mem.scale.prior,
      mem.overdispersion=mem.overdispersion,
      muproposal=muproposal, nuproposal=nuproposal, 
      beta0proposal=beta0proposal, beta1proposal=beta1proposal,
      memmuproposal=memmuproposal, memscaleproposal=memscaleproposal,
      visibility=visibility,
      Np=Np,
      samplesize=samplesize.parallel,burnin=burnin,interval=interval,
      burnintheta=burnintheta,
      burninbeta=burninbeta,
      priorsizedistribution=priorsizedistribution,
      mean.prior.size=mean.prior.size, sd.prior.size=sd.prior.size,
      mode.prior.sample.proportion=mode.prior.sample.proportion,
      median.prior.size=median.prior.size,
      mode.prior.size=mode.prior.size,
      quartiles.prior.size=quartiles.prior.size,
      effective.prior.df=effective.prior.df,
      alpha=alpha,
      seed=seed,
      supplied=supplied,
      num.recruits=nr[!remvalues],
      recruit.times=recruit.times[!remvalues],
      num.recruits2=nr2[!remvalues2],
      recruit.times2=recruit.times2[!remvalues2],
      max.coupons=max.coupons,
      maxbeta=maxbeta)
#
#   Process the results
#
    ### Snow returns a list of length parallel where each element is the return of each posfn
    ### Following loops combines the separate MCMC samples into 1 using rbind, creating a matrix
    Cret <- outlist[[1]]
    Nparallel <- length(outlist)
    Cret$samplesize <- samplesize
    for(i in (2 : Nparallel)){
     z <- outlist[[i]]
     Cret$sample <- rbind(Cret$sample,z$sample)
     if(visibility){
       Cret$vsample <- rbind(Cret$vsample,z$vsample)
       if(!is.null(s2)){
        Cret$vsample2 <- rbind(Cret$vsample2,z$vsample2)
       }
     }
     Cret$predictive.visibility.count<-Cret$predictive.visibility.count+z$predictive.visibility.count
     Cret$predictive.visibility<-Cret$predictive.visibility+z$predictive.visibility
     Cret$predictive.degree<-Cret$predictive.degree+z$predictive.degree
    }
    Cret$predictive.visibility.count<-Cret$predictive.visibility.count/Nparallel
    Cret$predictive.visibility<-Cret$predictive.visibility/Nparallel
    Cret$predictive.degree<-Cret$predictive.degree/Nparallel
    #
    degnames <- NULL
    if(Np>0){degnames <- c(degnames,paste("pdeg",1:Np,sep=""))}
#   colnamessample <- c("N","mu","sigma","visibility1","totalsize","beta0","beta1","mem.optimism","mem.scale","mem.visibility.mean")
#   colnamessample <- Cret$sample
#   if(length(degnames)>0){
#    colnamessample <- c(colnamessample, degnames)
#   }
#   if(visibilitydistribution=="cmp"){
#    colnamessample <- c(colnamessample, c("lambda","nu"))
#   }
#   colnames(Cret$sample) <- colnamessample
    m <- apply(Cret$sample,2,stats::median,na.rm=TRUE)
    Cret$sample[is.na(Cret$sample[,"mu"]),"mu"] <- m["mu"]
    Cret$sample[is.na(Cret$sample[,"sigma"]),"sigma"] <- m["sigma"]
#   Any NA and NaN are typically in pdeg and so should be 0.
    Cret$sample[is.na(Cret$sample)] <- 0
    Cret$sample[is.nan(Cret$sample)] <- 0
    
    ### Coda package which does MCMC diagnostics, requires certain attributes of MCMC sample
    endrun <- burnin+interval*(samplesize)
    attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    
#   ### Remove the padding from the last draws from the populations of visibilities
#   Nlastpos <- Cret$sample[nrow(Cret$sample),"N"]
#   Cret$pop<-Cret$pop[1:Nlastpop]

    ### compute modes of posterior samples,Maximum A Posterior (MAP) values N, mu, sigma, visibility1
    Cret$MAP <- apply(Cret$sample,2,mode.density)
    Cret$MAP["N"] <- mode.density(Cret$sample[,"N"],lbound=n,ubound=Cret$maxN)
    if(verbose){
     message("parallel samplesize = ", parallel," by ", samplesize.parallel,"\n",appendLF=FALSE)
    }
    
    ### stop cluster
    endparallel(cl,type=parallel.type)
  }
  Cret$N <- c(Cret$MAP["N"], 
              mean(Cret$sample[,"N"]),
              stats::median(Cret$sample[,"N"]),
              stats::quantile(Cret$sample[,"N"],c(0.025,0.975)))
  names(Cret$N) <- c("MAP","Mean AP","Median AP","P025","P975")
  #
  Cret$sample <- Cret$sample[,-match(c("visibility1","totalsize"), colnames(Cret$sample))]
  #
  if(verbose & Cret$predictive.visibility[length(Cret$predictive.visibility)] > 0.3){
   warning("There is a non-trivial proportion of the posterior mass on very high visibilities. This may indicate convergence problems in the MCMC.", call. = FALSE)
  }
  Cret$visibilitydistribution <- visibilitydistribution
  Cret$priorsizedistribution <- priorsizedistribution
  #
  if(missing(type.impute)){type.impute <- "median"}
  type.impute <- match.arg(type.impute,
                           c("median","distribution","mode","mean"))
  if(is.na(type.impute)) { # User typed an unrecognizable name
    stop(paste('You must specify a valid type.impute. The valid types are "distribution","mode","median", and "mean"'), call.=FALSE)
  }
  if(visibility){
    visibilities <- rep(0,length=nrow(rds.data))
    visibilities[!remvalues] <- switch(type.impute, 
               `distribution` = {
                 Cret$pop[1:Cret$n1]
               },
               `mode` = {
                 apply(Cret$vsample,2,function(x){a <- tabulate(x);mean(which(a==max(a,na.rm=TRUE)))})
               },
               `median` = {
                 apply(Cret$vsample,2,stats::median)
               },
               `mean` = {
                 apply(Cret$vsample,2,mean)
               }
              )
    # impute the missing values
    if(sum(remvalues) > 0){
     for(i in seq_along(recruit.times[remvalues])){
      mf <- (recruit.times[remvalues])[i] == recruit.times[!remvalues] & (nr[remvalues])[i] == nr[!remvalues]
      if(length(mf) > 0){
       mf <- matrix(Cret$vsample[,mf],ncol=1)
      }else{
       mf <- (recruit.times[remvalues])[i] == recruit.times[!remvalues]
       if(length(mf) > 0){
        mf <- matrix(Cret$vsample[,mf],ncol=1)
       }else{
        mf <- (nr[remvalues])[i] == nr[!remvalues]
        if(length(mf) > 0){
         mf <- matrix(Cret$vsample[,mf],ncol=1)
        }else{
         mf <- matrix(Cret$vsample,ncol=1)
        }
       }
      }
      visibilities[which(remvalues)[i]] <- switch(type.impute, 
                   `distribution` = {
                     mf[sample.int(n=1,size=length(mf))]
                   },
                   `mode` = {
                     apply(mf,2,function(x){a <- tabulate(x);mean(which(a==max(a,na.rm=TRUE)))})
                   },
                   `median` = {
                     apply(mf,2,stats::median)
                   },
                   `mean` = {
                     apply(mf,2,mean)
                   }
                 )
      }
    }
#     print(length(visibilities))
#     print(dim(Cret$vsample))
#     print(range(recruit.times.order))
#     print(range(recruit.times.order.notrem))
    Cret$visibilities <- visibilities[recruit.times.order]
    Cret$vsample <- Cret$vsample[,recruit.times.order.notrem]
    if(!is.null(s2)){
#     print(dim(Cret$vsample2))
#     print(range(recruit.times2.order))
#     print(range(recruit.times2.order.notrem))
      Cret$vsample2 <- Cret$vsample2[,recruit.times2.order.notrem]
    }
  }

  Cret$visibility <- visibility
# Cret$mean.prior.size <- mean.prior.size
  Cret$data <- rds.data
  Cret$data2 <- rds.data2
  ### return result
  class(Cret) <- "sspse"
  Cret
}

#' Warning message for posteriorsize fit failure
#' 
#' \code{\link{posteriorsize}} computes the posterior distribution of the
#' population size based on data collected by Respondent Driven Sampling.
#' This function returns the warning message if it fails. 
#' It enables packages that call \code{\link{posteriorsize}} to use
#' a consistent error message.
#' @return \code{\link{posize_warning}} returns a character string witn the warning message.
#' @seealso posteriorsize
#' @keywords models
#' @export posize_warning
posize_warning <- function(){
        "POSTERIOR SIZE CALCULATION FAILED" #added for posteriorsize dialog to hide error message unless needed
}
