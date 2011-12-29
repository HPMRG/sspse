dsizeprior<-function(n,
		  type=c("proportion","nbinom","pln","flat","continuous"),
		  mean.prior.size=NULL, sd.prior.size=NULL,
		  mode.prior.sample.proportion=0.5,
		  median.prior.sample.proportion=NULL,
		  median.prior.size=NULL,
		  median.mid.prior.size=NULL,
		  mode.prior.size=NULL,
		  effective.prior.df=1,
		  alpha=NULL,
                  maxN=NULL,
                  log=FALSE,
                  maxbeta=100,
                  maxNmax=200000,
                  verbose=TRUE){
  priorsizedistribution=match.arg(type)
  N <- NULL
  beta <- NULL
  lpriorm <- switch(priorsizedistribution,
    nbinom={
      if(is.null(sd.prior.size)){sd.prior.size <- mean.prior.size}
      if(is.null(maxN)){
        maxN <- min(maxNmax,ceiling(qnbinommu(p=0.995,
                    mu=mean.prior.size, sd=sd.prior.size)))
      }
      if(is.null(N)){
        maxN <- min(maxNmax,ceiling(qnbinommu(p=0.5,
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
     if(is.null(maxN)){maxN <- min(maxNmax,ceiling( n/(1-0.90^(1/beta)) ))}
     if(is.null(N)){N <- min(maxNmax,ceiling( n/(1-0.5^(1/beta)) ))}
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
     lfn <- function(x,beta,n,effective.prior.df,alpha){
#     (beta*n*((x-n)^(beta-1))/(x^(beta+1)))^effective.prior.df
#     a=effective.prior.df*(log(beta*n)+(beta-1)*log(x-n) - (beta+1)*log(x) )
      a=effective.prior.df*(log(n)+lgamma(alpha+beta)-lgamma(alpha)-lgamma(beta)+(beta-1)*log(x-n) - (alpha+beta)*log(x) )
      mina <- min(a,na.rm=TRUE)
#if(any(is.na(a))){browser()}
      a[is.na(a)] <- mina
      a[is.infinite(a)] <- mina
      a <- a - mina
     }
     if(!is.null(median.mid.prior.size)){
      if(median.prior.size < n){median.prior.size = n}
      if(median.prior.size < 750){effective.prior.df=max(effective.prior.df,3)}
      beta <- -log(2)/log(1-n/median.prior.size)
      if(is.null(alpha)) alpha=median.prior.size/(median.prior.size-n)
      if(is.null(maxN)){maxN <- min(maxNmax,ceiling( n/(1-0.90^(1/beta)) ))}
      fn <- function(beta,x,n,median.prior.size,effective.prior.df,alpha){
       lpriorm <- lfn(x+0.5,beta,n,effective.prior.df,alpha)
       priorm <- exp(lpriorm) / sum(exp(lpriorm) )
#   abs(median.prior.size - (x+0.5)[match(TRUE,cumsum(priorm) >= 0.5)] ) 
  print(c(beta,sum(x*priorm)/sum(priorm),
         0.5*( (x+0.5)[match(TRUE,cumsum(priorm) >= 0.5)] 
          +(x+0.5)[which.max(lpriorm)] )
        ))
       abs(median.prior.size - 
         0.5*( (x+0.5)[match(TRUE,cumsum(priorm) >= 0.5)] 
          +(x+0.5)[which.max(lpriorm)] )
          )
      }
      print(paste("median:",median.prior.size))
      while(pbeta(n/maxN,shape1=alpha,shape2=beta)<0.01){
       maxN = min(5000,qbeta(0.01,shape1=alpha,shape2=beta))
       x <- 0:(maxN-1) + n
       a = optimize(f=fn,interval=c(1,maxbeta),x,n,median.prior.size,
                    effective.prior.df,alpha,tol=0.01)
       beta <- a$minimum
      }
     }
     if(!is.null(median.prior.size)){
      if(median.prior.size < n){median.prior.size = n}
      if(median.prior.size < 750){effective.prior.df=max(effective.prior.df,3)}
      beta <- -log(2)/log(1-n/median.prior.size)
      if(is.null(alpha)) alpha=median.prior.size/(median.prior.size-n)
      if(is.null(maxN)){maxN <- min(maxNmax,ceiling( n/(1-0.90^(1/beta)) ))}
      fn <- function(beta,x,n,median.prior.size,effective.prior.df,alpha){
       lpriorm <- lfn(x+0.5,beta,n,effective.prior.df,alpha)
       priorm <- exp(lpriorm) / sum(exp(lpriorm) )
  print(c(beta,sum(x*priorm)/sum(priorm),
       (x+0.5)[match(TRUE,cumsum(priorm) >= 0.5)] ))
       abs(median.prior.size - 
        (x+0.5)[match(TRUE,cumsum(priorm) >= 0.5)] ) 
      }
  print(paste("median:",median.prior.size))
      maxN = ceiling(3*median.prior.size)
      x <- 0:(maxN-1-n) + n
      a = optimize(f=fn,interval=c(1,maxbeta),x,n,median.prior.size,
                   effective.prior.df,alpha,tol=0.01)
      beta <- a$minimum
      while(abs(pbeta(n/maxN,shape1=alpha,shape2=beta)-0.02)>0.001){
       maxN = ceiling(n/qbeta(0.02,shape1=alpha,shape2=beta))
       x <- 0:(maxN-1-n) + n
       a = optimize(f=fn,interval=c(1,maxbeta),x,n,median.prior.size,
                    effective.prior.df,alpha,tol=0.01)
       beta <- a$minimum
print(c(beta,maxN,pbeta(n/maxN,shape1=alpha,shape2=beta)))
      }
      maxN = min(maxNmax,maxN-n)
     }
     if(!is.null(mean.prior.size)){
      if(mean.prior.size < n){mean.prior.size = n}
      beta <- max(1.1,mean.prior.size/n - 1)
      if(is.null(alpha)) alpha=mean.prior.size/(mean.prior.size-n)
      beta <- 0.5*(1+alpha)
      if(is.null(maxN)){maxN <- min(maxNmax,ceiling( n/(1-0.90^(1/beta)) ))}
      fn <- function(beta,x,n,mean.prior.size,effective.prior.df,alpha){
#      lpriorm <- (1-n/(x+0.5))^beta - (1-n/(x-0.5))^beta
#      lpriorm[1] <- (1-n/(n+0.5))^beta
#      lpriorm <- lpriorm^effective.prior.df
       lpriorm <- lfn(x+0.5,beta,n,effective.prior.df,alpha)
#      lpriorm[is.infinite(lpriorm)] <- 0
#      lpriorm[is.na(lpriorm)] <- 0
#      lpriorm[is.na(lpriorm)] <- min(lpriorm,na.rm=TRUE)
       priorm <- exp(lpriorm)
#print(c(beta,sum(x*priorm)/sum(priorm),
#        abs(mean.prior.size - sum(x*priorm)/sum(priorm))))
       abs(mean.prior.size - sum(x*priorm)/sum(priorm,na.rm=TRUE))
      }
#     maxNfn <- function(maxN,n,mean.prior.size,effective.prior.df,alpha){
#      x <- 0:(maxN-1) + n
#      lpriorm <- lfn(x+0.5,beta,n,effective.prior.df,alpha)
#      priorm <- exp(lpriorm)
#      abs(0.99 - sum(priorm,na.rm=TRUE))
#     }
      maxN = ceiling(3*mean.prior.size)
      x <- 0:(maxN-1-n) + n
#     a = optimize(f=maxNfn,interval=c(1,maxNmax),n,mean.prior.size,
#                  effective.prior.df,alpha,tol=0.01)
      a = optimize(f=fn,interval=c(1,maxbeta),x,n,mean.prior.size,
                   effective.prior.df,alpha,tol=0.01)
      beta <- a$minimum
      while(abs(pbeta(n/maxN,shape1=alpha,shape2=beta)-0.02)>0.001){
       maxN = ceiling(n/qbeta(0.02,shape1=alpha,shape2=beta))
       x <- 0:(maxN-1-n) + n
       a = optimize(f=fn,interval=c(1,maxbeta),x,n,mean.prior.size,
                    effective.prior.df,alpha,tol=0.01)
       beta <- a$minimum
print(c(beta,maxN,pbeta(n/maxN,shape1=alpha,shape2=beta)))
      }
      maxN = min(maxNmax,maxN)
     }
     if(!is.null(mode.prior.size)){
      if(mode.prior.size < n){mode.prior.size = n}
      beta <- 2*mode.prior.size/n - 1
      if(is.null(alpha)) alpha=mode.prior.size/(mode.prior.size-n)
      if(is.null(maxN)){maxN <- min(maxNmax,ceiling( n/(1-0.90^(1/beta)) ))}
      fn <- function(beta,x,n,mode.prior.size,effective.prior.df,alpha){
       lpriorm <- lfn(x+0.5,beta,n,effective.prior.df,alpha)
       priorm <- exp(lpriorm)
  print(c(beta,sum(x*priorm)/sum(priorm),
          abs(mode.prior.size - (x+0.5)[which.max(lpriorm)]) ))
       abs(mode.prior.size - (x+0.5)[which.max(lpriorm)])
      }
      maxN = ceiling(3*mode.prior.size)
      x <- 0:(maxN-1-n) + n
      a = optimize(f=fn,interval=c(1,maxbeta),x,n,mode.prior.size,
                   effective.prior.df,alpha,tol=0.01)
      beta <- a$minimum
      while(abs(pbeta(n/maxN,shape1=alpha,shape2=beta)-0.02)>0.001){
       maxN = ceiling(n/qbeta(0.02,shape1=alpha,shape2=beta))
       x <- 0:(maxN-1-n) + n
       a = optimize(f=fn,interval=c(1,maxbeta),x,n,mode.prior.size,
                    effective.prior.df,alpha,tol=0.01)
       beta <- a$minimum
print(c(beta,maxN,pbeta(n/maxN,shape1=alpha,shape2=beta)))
      }
      maxN = min(maxNmax,maxN)
     }
#    median.prior.size <- n/(1-0.5^(1/beta))
#    mode.prior.size <- n*(beta+1)/2
#    mean.prior.size <- n*(beta+1)
#    mode.prior.sample.proportion <- 2/(beta+1)
     if(is.null(maxN)){maxN <- min(maxNmax,ceiling( n/(1-0.90^(1/beta)) ))}
     if(is.null(N)){N <- min(maxNmax,ceiling( n/(1-0.5^(1/beta)) ))}
#    x <- 0:maxN
#    a <- 0.5*(n/(n+x) + n/(n+x+1))
#    priorm <- (1-a)^beta
#    lpriorm <- diff(priorm)
#    lpriorm[1] <- priorm[1]
     x <- 0:(maxN-1) + n
#    lpriorm <- (1-n/(x+0.5))^beta - (1-n/(x-0.5))^beta
#    lpriorm[1] <- (1-n/(n+0.5))^beta
#    lpriorm <- ((beta*n)^effective.prior.df)*( 
#       (((x+0.5-n)^(beta-1))/((x+0.5)^(beta+1)))^effective.prior.df -
#       (((x-0.5-n)^(beta-1))/((x-0.5)^(beta+1)))^effective.prior.df ) 
#    lpriorm[1] <- ((beta*n)*(0.5)^(beta-1)/(n+0.5)^(beta+1))^effective.prior.df
#      lpriorm <-
#      lfn(x+0.5,beta,n,effective.prior.df)-lfn(x-0.5,beta,n,effective.prior.df,alpha)
#      lpriorm[1] <- lfn(n+0.5,beta,n,effective.prior.df,alpha)
     lpriorm <- lfn(x+0.5,beta,n,effective.prior.df,alpha)
     priorm <- exp(lpriorm)
     lpriorm <- priorm / sum(priorm,na.rm=TRUE)
     if(log){
      lpriorm <- log(lpriorm)
     }
     x <- 0:(maxN-1)
     lpriorm
     }
    )
#   if(effective.prior.df!=1){
     if(log){
#    lpriorm <- exp(lpriorm*effective.prior.df-max(lpriorm*effective.prior.df,na.rm=TRUE))
      priorm <- exp(lpriorm)
#     mean.prior.size <- sum((x+n+0.5)*priorm)/sum(priorm)
#     median.prior.size <- (x+n+0.5)[match(TRUE,cumsum(priorm) >= sum(priorm)*0.5)]
#     mode.prior.size <- (x+n+0.5)[which.max(priorm)]
     }else{
      priorm <- lpriorm
     }
     mean.prior.size <- sum((x+n+0.5)*priorm)
     median.prior.size <- (x+n+0.5)[match(TRUE,cumsum(priorm) >= 0.5)]
     mode.prior.size <- (x+n+0.5)[which.max(priorm)]
#   }
    if(is.null(N)){N <- mean(x)}
    if(verbose){cat(paste("The maximum prior population size is",n+maxN,"\n"))}
    list(x=n+x,lprior=lpriorm,N=N,maxN=maxN,
         median.prior.size=median.prior.size,
         mean.prior.size=mean.prior.size,
         mode.prior.size=mode.prior.size,
	 mode.prior.sample.proportion=mode.prior.sample.proportion,
	 median.prior.sample.proportion=median.prior.sample.proportion,
	 alpha=alpha,beta=beta,
	 effective.prior.df=effective.prior.df,
         type=type)
}
