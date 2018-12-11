#
# Calculate the Measurement error model using Bayesian Inference
#
#
#' Estimates each person's personal visibility based on their self-reported degree and the 
#'  number of their (direct) recruits. It uses the time the person was recruited as a factor in 
#'  determining the number of recruits they produce.
#' @param rds.data An rds.data.frame
#' @param max.coupons The number of recruitment coupons distributed to each 
#' 		enrolled subject (i.e. the maximum number of recruitees for any subject).
#'              By default it is taken by the attribute or data, else the maximum recorded number of coupons.
#' @param type.impute The type of imputation based on the conditional distribution. 
#'     It can be of type \code{distribution},\code{mode},\code{median}, or \code{mean} 
#'     with the first , the default, being a random draw from the conditional distribution.
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
#' @param unit.scale numeric; If not \code{NULL} it sets the numeric value of the scale parameter
#' of the distribution of the unit sizes.
#' For the negative binomial, it is the multiplier on the variance of the negative binomial 
#' compared to a Poisson (via the Poisson-Gamma mixture representation). Sometimes the scale is 
#' unnaturally large (e.g. 40) so this give the option of fixing it (rather than using
#' the MLE of it). The model is fit with the parameter fixed at this passed value.
#' It can be of \code{nbinom}, meaning a negative binomial. 
#' In this case, \code{unit.scale} is the multiplier 
#' on the variance of the negative binomial compared to a Poisson of the same mean.
#' The alternative is \code{cmp}, meaning a Conway-Maxwell-Poisson distribution.
#' In this case, \code{unit.scale}
#' is the scale parameter compared to a Poisson of the same mean (values less than one mean 
#' under-dispersed and values over one mean over-dispersed). The default is \code{cmp}.
#' @param reflect.time logical; If \code{FALSE} then the \code{recruit.time} is the time before the 
#' end of the study (instead of the time since the survey started or chronological time).
#' @param K count; the maximum visibility for an individual. This is usually
#' calculated as \code{round(stats::quantile(s,0.80))}. It applies to network sizes and (latent) visibilities.
#' If logical and FALSE then the K is unbounded but set to compute the visibilities.
#' @param parallel count; the number of parallel processes to run for the
#' Monte-Carlo sample.  This uses MPI or PSOCK. The default is 1, that is not to
#' use parallel processing.
#' @param parallel.type The type of parallel processing to use. The options are
#' "PSOCK" or "MPI". This requires the corresponding type to be installed.
#' The default is "PSOCK".
#' @param interval count; the number of proposals between sampled statistics.
#' @param burnin count; the number of proposals before any MCMC sampling is
#' done. It typically is set to a fairly large number.
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
#' McLaughlin, K.R., M.S. Handcock, and L.G. Johnston, 2015. 
#' Inference for the visibility distribution for respondent-driven sampling. 
#' In JSM Proceedings. Alexandria, VA: American Statistical Association. 2259-2267.
impute.visibility <-function(rds.data,max.coupons=NULL,
                             type.impute = c("median","distribution","mode","mean"),
                             recruit.time=NULL,include.tree=FALSE, unit.scale=NULL, 
                             reflect.time=TRUE,
                             K=FALSE,
                             parallel=1, parallel.type="PSOCK",
                             interval=10,
                             burnin=5000,
                             verbose=TRUE){
  if(!is(rds.data,"rds.data.frame"))
    stop("rds.data must be of type rds.data.frame")   
  
  if(missing(type.impute)){type.impute <- "median"}

  fit <- posteriorsize(s=rds.data,
                  K=K,
                  df.mean.prior.visibility=1000,
                  visibility=TRUE,
		  type.impute = type.impute,
                  parallel=parallel, parallel.type=parallel.type,
                  recruit.time=recruit.time,
		  include.tree=include.tree, unit.scale=unit.scale, 
                  optimism = TRUE,
                  reflect.time=reflect.time,
                  verbose=TRUE)

  if(verbose){
#   print(summary(fit))
    print(prettyNum(apply(fit$sample,2,median)[-1],digits=2))
  }
  
  return(round(fit$visibilities))
}
