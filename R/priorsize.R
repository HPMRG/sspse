priorsize<-function(s,
		  median.prior.size=NULL,
                  maxN=NULL,
                  K=NULL,
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
                  Np=0,
                  nk=NULL,
                  n=length(s),
                  seed=NULL, 
                  supplied=list(maxN=maxN),
                  verbose=TRUE){
#
  visibilitydistribution=match.arg(visibilitydistribution)
  posfn <- switch(visibilitydistribution,
                  nbinom=posnbinom,
                  pln=pospln,
		  cmp=poscmp,
		  poscmp)
  priorsizedistribution=match.arg(priorsizedistribution)
  if(priorsizedistribution=="nbinom" && missing(mean.prior.size)){
    stop("You need to specify 'mean.prior.size', and possibly 'sd.prior.size' if you use the 'nbinom' prior.") 
  }
  if(is.null(K)){
    K=round(stats::quantile(s,0.80))
    degs <- s
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
    if(sd.pd > sqrt(max.sd.prior.visibility*mean.pd)){
     sd.pd <- min(sqrt(max.sd.prior.visibility*mean.pd), sd.pd)
    }
    xv <- ds
#   xp <- weights*ds
    xp <- weights
    xp <- length(xp)*xp/sum(xp)
    fit <- cmpmle(xv,xp,cutoff=1,cutabove=K-1,guess=c(mean.pd, sd.pd))
    y=dcmp.natural(v=fit,x=(0:max(s)))
    K=(0:max(s))[which.max(cumsum(y)>0.95)]
  }
  cat(sprintf("The size cap is K = %d.\n",K))
  if(is.null(mean.prior.visibility)){
    degs <- s
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
    xp <- length(xp)*xp/sum(xp)
    fit <- cmpmle(xv,xp,cutoff=1,cutabove=K-1,
            guess=c(mean.prior.visibility,sd.prior.visibility))
    fit <- cmp.to.mu.sd(fit,max.mu=5*mean.prior.visibility)
    mean.prior.visibility = fit[1]
    sd.prior.visibility = fit[2]
  }
  if(verbose){
    cat(sprintf("The mean of the prior distribution for visibility is %f.\n",mean.prior.visibility))
    cat(sprintf("The s.d. of the prior distribution for visibility is %f.\n",sd.prior.visibility))
  }
  if(sd.prior.visibility > sqrt(max.sd.prior.visibility*mean.prior.visibility)){
    sd.prior.visibility <- min(sqrt(max.sd.prior.visibility*mean.prior.visibility), sd.prior.visibility)
    cat(sprintf("The suggested s.d. of the prior distribution for visibility is too large and has been reduced to the more reasonable %f.\n",sd.prior.visibility))
  }
  ### are we running the job in parallel (parallel > 1), if not just 
  #   call the visibility specific function
    if(!is.null(seed))  set.seed(as.integer(seed))
    #
    # Cap the maximum visibility to K
    #
    s[s>K] <- K
    if(is.null(nk)){nk=tabulate(s,nbins=K)}
    #
    # Transform observed mean parametrization to log-normal
    # parametrization
    #
    out <- cmp.to.natural(mu=mean.prior.visibility, sigma=sd.prior.visibility)
    mu <- log(out$lambda)
    sigma <- out$nu
    #
    dimsample <- 5+Np
    #
    priorsizedistribution=match.arg(priorsizedistribution)
    prior <- dsizeprior(n=n,
		  type=priorsizedistribution,
		  sd.prior.size=sd.prior.size,
		  mode.prior.sample.proportion=mode.prior.sample.proportion,
		  median.prior.size=median.prior.size,
		  mode.prior.size=mode.prior.size,
		  mean.prior.size=mean.prior.size,
		  quartiles.prior.size=quartiles.prior.size,
                  effective.prior.df=effective.prior.df,
                  alpha=alpha,
                  maxN=maxN,
                  log=FALSE,
                  supplied=supplied,
                  verbose=verbose)
   x <- seq_along(prior) + n - 1

   M <- 1000
   mx <- rep(0,M)
   sdx <- rep(0,M)
   nk <- rep(0,K)
   for(m in 1:M){
    N <- sample(prior$x,size=1,prob=prior$lprior)
    if(N > n){
     sigma <- sd.prior.visibility*sqrt(1/(stats::rchisq(n=1, df = df.sd.prior.visibility)))
     mu <- stats::rnorm(n=1,mean=mean.prior.visibility, sd=sigma / df.mean.prior.visibility)
     pcmp <- cmp.to.natural(mu=mu,sigma=sigma)
     pcmp <- .C("rcmp",
                x=integer(N-n),
                lambda=as.double(pcmp$lambda),
                nu=as.double(pcmp$nu),
                n=as.integer(N-n),
                K=as.integer(K),
                err=as.double(0.00000001),
                PACKAGE="sspse")$x
     nk=nk+( tabulate(c(s,pcmp),nbins=K) / (sum(nk > 0) ) )
     mx[m] <- mean(c(s,pcmp))
     sdx[m] <- sqrt(stats::var(c(s,pcmp)))
    }
   }
   list(mu=mean(mx),sd=mean(sdx),pmf=nk / M)
}
