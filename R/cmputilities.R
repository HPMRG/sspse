dsizeprior<-function(n,
		  type=c("proportion","nbinom","pln","flat"),
		  mean.prior.size=NULL, sd.prior.size=NULL,
		  mode.prior.sample.proportion=0.5,
		  median.prior.sample.proportion=NULL,
		  median.prior.size=NULL,
		  mode.prior.size=NULL,
                  maxN=NULL,
                  log=FALSE,
                  verbose=TRUE){
  priorsizedistribution=match.arg(type)
  if(priorsizedistribution=="nbinom" && is.null(mean.prior.size)){
    mean.prior.size<-N
  }
  N <- NULL
  beta <- NULL
  lpriorm <- switch(priorsizedistribution,
    nbinom={
      if(is.null(sd.prior.size)){sd.prior.size <- mean.prior.size}
      if(is.null(maxN)){
        maxN <- min(50000,ceiling(qnbinommu(p=0.995,
                    mu=mean.prior.size, sd=sd.prior.size)))
      }
      if(is.null(N)){
        maxN <- min(50000,ceiling(qnbinommu(p=0.5,
                    mu=mean.prior.size, sd=sd.prior.size)))
      }
      lpriorm <- dnbinommu(x=n+(1:maxN)-1,
                           mu=mean.prior.size, sd=sd.prior.size,
                           log=log)
      lpriorm
     },
    flat={
      if(is.null(maxN)){
        maxN <- 10*length(s)
      }
      if(is.null(N)){
        N <- 0.5*maxN
      }
      if(log){
       lpriorm <- rep(0,maxN)
      }else{
       lpriorm <- rep(1/maxN,maxN)
      }
      lpriorm
     },
    proportion={
     if(!is.null(mode.prior.sample.proportion)){
      beta <- 2/mode.prior.sample.proportion - 1
     }
     if(!is.null(median.prior.sample.proportion)){
      beta <- -log(2)/log(1-median.prior.sample.proportion)
     }
     if(!is.null(median.prior.size)){
      beta <- -log(2)/log(1-n/median.prior.size)
     }
     if(!is.null(mean.prior.size)){
      beta <- (1-mean.prior.size)/mean.prior.size
     }
     if(!is.null(mode.prior.size)){
      beta <- 2*mode.prior.size/n - 1
     }
     median.prior.size <- n/(1-0.5^(1/beta))
     mode.prior.size <- n*(beta+1)/2
     mode.prior.sample.proportion <- 2/(beta+1)
     if(is.null(maxN)){maxN <- min(50000,ceiling( n/(1-0.90^(1/beta)) ))}
     if(is.null(N)){N <- min(50000,ceiling( n/(1-0.5^(1/beta)) ))}
     x <- (1:maxN)
     lpriorm <- log(beta*n)+(beta-1)*log(x)-(beta+1)*log(x+n)
     if(!log){
      lpriorm <- exp(lpriorm)
     }
     lpriorm
     }
    )
    if(verbose){cat(paste("The maximum prior range is",n+maxN,"\n"))}
    list(x=n+x,lprior=lpriorm,N=N,maxN=maxN,
         median.prior.size=median.prior.size,
         mean.prior.size=mean.prior.size,
         mode.prior.size=mode.prior.size,
	 mode.prior.sample.proportion=mode.prior.sample.proportion,
	 median.prior.sample.proportion=median.prior.sample.proportion,
	 beta=beta)
}

rcmp.lambda <- function(n, lambda, nu, err=0.00001, K=100){
  # Perform argument checking
  if (lambda < 0 || nu < 0)
	stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
  if (n < 0 || n != floor(n))
	stop("Invalid number of draws");
  out <- .C("rcmp",
            x=integer(n),
            lambda=as.double(lambda),
            nu=as.double(nu),
            n=as.integer(n),
            K=as.integer(K),
            err=as.double(err),
            PACKAGE="size")
   return(out$x)
}

rcmp <- function(n, mu, sig, err=0.00001, K=100){
  out <- cmp.natural(mu,sig)
  # Perform argument checking
  if (out$lambda < 0 || out$nu < 0)
	stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
  if (n < 0 || n != floor(n))
	stop("Invalid number of draws");
  out <- .C("rcmp",
            x=integer(n),
            lambda=as.double(out$lambda),
            nu=as.double(out$nu),
            n=as.integer(n),
            K=as.integer(K),
            err=as.double(err),
            PACKAGE="size")
   return(out$x)
}

dcmp <- function(x, lambda, nu, err=0.00001, log=FALSE){
  # Perform argument checking
  if (lambda < 0 || nu < 0)
	stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
  if (any(x < 0 || x != floor(x)))
		return (0);
  out <- .C("dcmp",
            x=as.integer(x),
            lambda=as.double(lambda),
            nu=as.double(nu),
            n=as.integer(length(x)),
            err=as.double(err),
            give_log=as.integer(log),
            val=double(length(x)),
            PACKAGE="size")
   return(out$val)
}

dcmp.mu <- function(x, mu, sig, err=0.00001, log=FALSE){
  # Perform argument checking
  if (mu < 0 || sig < 0)
	stop("Invalid arguments, only defined for mu >= 0, sd >= 0");
  if (any(x < 0 || x != floor(x)))
		return (0);
  out <- cmp.natural(mu,sig)
  out <- .C("dcmp",
            x=as.integer(x),
            lambda=as.double(out$lambda),
            nu=as.double(out$nu),
            n=as.integer(length(x)),
            err=as.double(err),
            give_log=as.integer(log),
            val=double(length(x)),
            PACKAGE="size")
   return(out$val)
}
cmp.compute.z = function(lambda, nu, error = 0.001, log=FALSE)
{
    if (missing(lambda)) {
      stop("The mean of the Negative Binomial must be specified.")
    }
    if (missing(nu)) {
      stop("The standard deviation of the Negative Binomial must be specified.")
    }
#  return(out$val)
   out <- .C("vzcmp",
            lambda=as.double(lambda),
            nu=as.double(nu),
            err=as.double(error),
            give_log=as.integer(log),
            out=as.double(error),
            PACKAGE="size")
   return(out$out)
}

dcomsize = function(x, lambda, nu, z = NULL)
{
	# Perform argument checking
	if (lambda < 0 || nu < 0)
		stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
	if (x < 0 || x != floor(x))
		return (0);
	if (is.null(z) || z <= 0)
		z = cmp.compute.z(lambda, nu);
	
	# Return pmf
	return ((lambda ^ x) * ((factorial(x)) ^ -nu) / z);
}
