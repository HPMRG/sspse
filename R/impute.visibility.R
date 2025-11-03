#
# Estimate each person's visibility based on a measurement error model and Bayesian inference
#
#
#' Estimates each person's personal visibility based on their self-reported degree and the 
#' number of their (direct) recruits. It uses the time the person was recruited as a factor in 
#' determining the number of recruits they produce.
#' @param rds.data An rds.data.frame
#' @param max.coupons The number of recruitment coupons distributed to each 
#' enrolled subject (i.e. the maximum number of recruitees for any subject).
#' By default it is taken by the attribute or data, else the maximum recorded number of coupons.
#' @param type.impute The type of imputation based on the conditional distribution. 
#' It can be of type \code{distribution},\code{mode},\code{median}, or \code{mean} 
#' with the first , the default, being a random draw from the conditional distribution.
#' @param recruit.time vector; An optional value for the data/time that the person was interviewed.
#' It needs to resolve as a numeric vector with number of elements the number
#' of rows of the data with non-missing values of the network variable. If it
#' is a character name of a variable in the data then that variable is used.
#' If it is NULL then the sequence number of the recruit in the data is used.
#' If it is NA then the recruitment is not used in the model.
#' Otherwise, the recruitment time is used in the model to better predict the visibility of the person.
#' @param include.tree logical; If \code{TRUE}, 
#' augment the reported network size by the number of recruits and one for the recruiter (if any).
#' This reflects a more accurate value for the visibility, but is not the self-reported degree.
#' In particular, it typically produces a positive visibility (compared to a possibility zero self-reported degree). 
#' @param reflect.time logical; If \code{FALSE} then the \code{recruit.time} is the time before the 
#' end of the study (instead of the time since the survey started or chronological time).
#' @param parallel count; the number of parallel processes to run for the
#' Monte-Carlo sample.  This uses MPI or PSOCK. The default is 1, that is not to
#' use parallel processing.
#' @param parallel.type The type of parallel processing to use. The options are
#' "PSOCK" or "MPI". This requires the corresponding type to be installed.
#' The default is "PSOCK".
#' @param samplesize count; the number of Monte-Carlo samples to draw to
#' compute the posterior. This is the number returned by the
#' Metropolis-Hastings algorithm. The default is 1000.
#' @param interval count; the number of proposals between sampled statistics.
#' @param warmup count; the number of proposals before any MCMC sampling is
#' done. It typically is set to a fairly large number.
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
#' @param return.posterior.sample.visibilities logical; If TRUE then return a
#' matrix of dimension \code{samplesize} by \code{n} of 
#' posterior draws from the visibility distribution for those in the survey.
#' The sample for the \code{i}th person is the \code{i}th column.
#' The default is FALSE so that the vector of imputes defined by \code{type.impute} is returned.
#' @param verbose logical; if this is \code{TRUE}, the program will print out additional
#  information about the fitting process.
#' @export impute.visibility
#' @examples
#' \dontrun{
#' data(fauxmadrona)
#' # The next line fits the model for the self-reported personal
#' # network sizes and imputes the personal network sizes 
#' # It may take up to 60 seconds.
#' visibility <- impute.visibility(fauxmadrona)
#' # frequency of estimated personal visibility
#' table(visibility)
#' }
#' @references 
#' McLaughlin, Katherine R.; Johnston, Lisa G.; Jakupi, Xhevat; Gexha-Bunjaku, Dafina; Deva, Edona and Handcock, Mark S. (2023)
#' Modeling the Visibility Distribution for Respondent-Driven Sampling with Application to Population Size Estimation,
#' \emph{Annals of Applied Statistics}, \doi{10.1093/jrsssa/qnad031} 
impute.visibility <-function(rds.data,max.coupons=NULL,
                             type.impute = c("median","distribution","mode","mean"),
                             recruit.time=NULL,include.tree=FALSE,
                             reflect.time=FALSE,
                             parallel=1, parallel.type="PSOCK",
                             interval=10,
                             warmup=5000,
                             samplesize=1000,
                             mem.optimism.prior=NULL, df.mem.optimism.prior=5,
                             mem.scale.prior=2, df.mem.scale.prior=10,
                             mem.overdispersion=15,
                             return.posterior.sample.visibilities=FALSE,
                             verbose=FALSE){
  if(!is(rds.data,"rds.data.frame"))
    stop("rds.data must be of type rds.data.frame")   
  
  if(missing(type.impute)){type.impute <- "median"}

    posfn <- posteriorsize
    fit <- posfn(s=rds.data,
               K=FALSE,
               visibility=TRUE,
               type.impute = type.impute,
               parallel=parallel, parallel.type=parallel.type,
               samplesize=samplesize, interval=interval, warmup=warmup,
               recruit.time=recruit.time, include.tree=include.tree,
               optimism = TRUE,
               reflect.time=reflect.time,
               mem.optimism.prior=mem.optimism.prior, df.mem.optimism.prior=df.mem.optimism.prior,
               mem.scale.prior=mem.scale.prior, df.mem.scale.prior=df.mem.scale.prior,
               mem.overdispersion=mem.overdispersion,
               verbose=verbose)

  if(verbose){
    print(prettyNum(apply(fit$sample,2,median)[-1],digits=2))
  }
  
  if(return.posterior.sample.visibilities){
   return(round(fit$vsample))
  }else{
   return(round(fit$visibilities))
  }
}
