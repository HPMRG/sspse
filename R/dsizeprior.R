dsizeprior<-function(n,
		  type=c("proportion","nbinom","pln","flat","continuous"),
		  mean.prior.size=NULL, sd.prior.size=NULL,
		  mode.prior.sample.proportion=0.5,
		  median.prior.sample.proportion=NULL,
		  median.prior.size=NULL,
		  mode.prior.size=NULL,
		  effective.prior.df=1,
                  maxN=NULL,
                  log=FALSE,
                  verbose=TRUE){
  priorsizedistribution=match.arg(type)
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
      x <- 0:(maxN-1)
      if(is.null(median.prior.size)) median.prior.size <- maxN/2
      lpriorm
     },
    flat={
      if(is.null(maxN)){
        maxN <- 10*n
      }
      if(is.null(N)){
        N <- 0.5*maxN
      }
      if(is.null(mode.prior.size)){
        mode.prior.size <- 0
      }
      if(log){
       lpriorm <- rep(0,maxN)
      }else{
       lpriorm <- rep(1/maxN,maxN)
       if(mode.prior.size > n){
        lpriorm[1:(mode.prior.size-n)] <- 0
       }
      }
      x <- 0:(maxN-1)
      if(is.null(median.prior.size)) median.prior.size <- maxN/2
      lpriorm
     },
    continuous={
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
      beta <- mean.prior.size/n - 1
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
      beta <- mean.prior.size/n - 1
     }
     if(!is.null(mode.prior.size)){
      beta <- 2*mode.prior.size/n - 1
     }
     median.prior.size <- n/(1-0.5^(1/beta))
     mode.prior.size <- n*(beta+1)/2
     mode.prior.sample.proportion <- 2/(beta+1)
     if(is.null(maxN)){maxN <- min(50000,ceiling( n/(1-0.90^(1/beta)) ))}
     if(is.null(N)){N <- min(50000,ceiling( n/(1-0.5^(1/beta)) ))}
#    x <- 0:maxN
#    a <- 0.5*(n/(n+x) + n/(n+x+1))
#    priorm <- (1-a)^beta
#    lpriorm <- diff(priorm)
#    lpriorm[1] <- priorm[1]
     x <- 0:(maxN-1) + n
     lpriorm <- (1-n/(x+0.5))^beta - (1-n/(x-0.5))^beta
     lpriorm[1] <- (1-n/(n+0.5))^beta
     if(log){
      lpriorm <- log(lpriorm)
     }
     x <- 0:(maxN-1)
     lpriorm
     }
    )
    if(effective.prior.df!=1){
     if(log){
      lpriorm <- exp(lpriorm*effective.prior.df)
      lpriorm <- lpriorm - max(lpriorm)
      lpriorm <- lpriorm/sum(lpriorm)
      lpriorm <- log(lpriorm)
      lpriorm[is.infinite(lpriorm)] <- NA
      lpriorm[is.na(lpriorm)] <- min(lpriorm,na.rm=TRUE)
     }else{
      lpriorm <- lpriorm^effective.prior.df
      lpriorm <- lpriorm/sum(lpriorm)
      lpriorm[is.infinite(lpriorm)] <- NA
      lpriorm[is.na(lpriorm)] <- min(lpriorm,na.rm=TRUE)
     }
    }
    if(is.null(N)){N <- mean(x)}
    if(verbose){cat(paste("The maximum prior population size is",n+maxN,"\n"))}
    list(x=n+x,lprior=lpriorm,N=N,maxN=maxN,
         median.prior.size=median.prior.size,
         mean.prior.size=mean.prior.size,
         mode.prior.size=mode.prior.size,
	 mode.prior.sample.proportion=mode.prior.sample.proportion,
	 median.prior.sample.proportion=median.prior.sample.proportion,
	 beta=beta,
	 effective.prior.df=effective.prior.df,
         type=type)
}
