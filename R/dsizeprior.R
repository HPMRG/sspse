dsizeprior<-function(n,
		  type=c("proportion","nbinom","pln","flat","continuous"),
		  mean.prior.size=NULL, sd.prior.size=NULL,
		  mode.prior.sample.proportion=NULL,
		  median.prior.sample.proportion=0.5,
		  median.prior.size=NULL,
		  median.mid.prior.size=NULL,
		  mode.prior.size=NULL,
		  quartiles.prior.size=NULL,
		  effective.prior.df=1,
		  alpha=NULL,
		  beta=NULL,
                  maxN=NULL,
                  log=FALSE,
                  maxbeta=100,
                  maxNmax=200000,
                  verbose=TRUE){
  priorsizedistribution=match.arg(type)
  N <- NULL
  lfn <- function(x,beta,n,effective.prior.df,alpha){
   a=effective.prior.df*(log(n)+lgamma(alpha+beta)-lgamma(alpha)-lgamma(beta)+(beta-1)*log(x-n) - (alpha+beta)*log(x) )
   mina <- min(a,na.rm=TRUE)
   a[is.na(a)] <- mina
   a[is.infinite(a)] <- mina
   a <- a - mina
   a
  }
  dfn <- function(alpha,beta,x,n,effective.prior.df){
   lpriorm <- lfn(x+0.5,beta,n,effective.prior.df,alpha)
   priorm <- exp(lpriorm)
   priorm/sum(priorm,na.rm=TRUE)
  }
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
      lpriorm <- dnbinommu(x=n:maxN,
                           mu=mean.prior.size, sd=sd.prior.size,
                           log=log)
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
       lpriorm <- rep(0,maxN-n+1)
      }else{
       lpriorm <- rep(1/maxN,maxN+1)
       if(mode.prior.size > n){
        lpriorm[1:(mode.prior.size-n)] <- 0
       }
      }
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
     x <- (n:maxN)
     lpriorm <- log(beta*n)+(beta-1)*log(x-n+1)-(beta+1)*log(x+1)
     if(!log){
      lpriorm <- exp(lpriorm)
     }
     lpriorm
     },
    proportion={
     if(is.null(alpha) | is.null(beta)){
     if(!is.null(mode.prior.sample.proportion)){
      beta <- 2/mode.prior.sample.proportion - 1
     }
     if(!is.null(median.prior.sample.proportion)){
      beta <- -log(2)/log(1-median.prior.sample.proportion)
     }
     if(!is.null(median.mid.prior.size)){
      if(median.prior.size < n){median.prior.size = n}
      if(median.prior.size < 750){effective.prior.df=max(effective.prior.df,3)}
      beta <- -log(2)/log(1-n/median.prior.size)
      if(is.null(alpha)) alpha=median.prior.size/(median.prior.size-n)
      if(is.null(maxN)){maxN <- min(maxNmax,ceiling( n/(1-0.90^(1/beta)) ))}
      fn <- function(beta,x,n,median.prior.size,effective.prior.df,alpha){
       priorm <- dfn(alpha,beta,x,n,effective.prior.df)
# print(c(beta,sum(x*priorm)/sum(priorm),
#        0.5*( (x+0.5)[match(TRUE,cumsum(priorm) >= 0.5)] 
#         +(x+0.5)[which.max(priorm)] )
#       ))
       abs(median.prior.size - 
         0.5*( (x+0.5)[match(TRUE,cumsum(priorm) >= 0.5)] 
          +(x+0.5)[which.max(priorm)] )
          )
      }
#     print(paste("median:",median.prior.size))
      x <- n:maxN
      while( {
       p=dfn(alpha,beta,x,n,effective.prior.df);
#      print( c(sum(p[x>1000]), sum(p[x>4000])) );
       abs(p[length(p)]/max(p) - 0.01)>0.005}){
        maxN <- round(maxN*(c(0.9,1.1)[(p[length(p)]/max(p) > 0.01)+1]))
        x <- n:maxN
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
       priorm <- dfn(alpha,beta,x,n,effective.prior.df)
# print(c(beta,sum(x*priorm)/sum(priorm),
#      (x+0.5)[match(TRUE,cumsum(priorm) >= 0.5)] ))
       abs(median.prior.size - 
        (x+0.5)[match(TRUE,cumsum(priorm) >= 0.5)] ) 
      }
# print(paste("median:",median.prior.size))
      maxN = ceiling(3*median.prior.size)
      x <- n:maxN
      a = optimize(f=fn,interval=c(1,maxbeta),x,n,median.prior.size,
                   effective.prior.df,alpha,tol=0.01)
      beta <- a$minimum
      while( {
       p=dfn(alpha,beta,x,n,effective.prior.df);
#      print( c(sum(p[x>1000]), sum(p[x>4000])) );
       abs(p[length(p)]/max(p) - 0.01)>0.005}){
        maxN <- round(maxN*(c(0.9,1.1)[(p[length(p)]/max(p) > 0.01)+1]))
        x <- n:maxN
        a = optimize(f=fn,interval=c(1,maxbeta),x,n,median.prior.size,
                     effective.prior.df,alpha,tol=0.01)
        beta <- a$minimum
#print(c(beta,maxN,pbeta(n/maxN,shape1=alpha,shape2=beta)))
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
       priorm <- dfn(alpha,beta,x,n,effective.prior.df)
       abs(mean.prior.size - sum(x*priorm)/sum(priorm,na.rm=TRUE))
      }
      maxN = ceiling(3*mean.prior.size)
      x <- n:maxN
      a = optimize(f=fn,interval=c(1,maxbeta),x,n,mean.prior.size,
                   effective.prior.df,alpha,tol=0.01)
      beta <- a$minimum
      while( {
       p=dfn(alpha,beta,x,n,effective.prior.df);
#      print( c(sum(p[x>1000]), sum(p[x>4000])) );
       abs(p[length(p)]/max(p) - 0.01)>0.005}){
        maxN <- round(maxN*(c(0.9,1.1)[(p[length(p)]/max(p) > 0.01)+1]))
        x <- n:maxN
        a = optimize(f=fn,interval=c(1,maxbeta),x,n,mean.prior.size,
                    effective.prior.df,alpha,tol=0.01)
        beta <- a$minimum
#print(c(beta,maxN,pbeta(n/maxN,shape1=alpha,shape2=beta)))
      }
      maxN = min(maxNmax,maxN)
     }
     if(!is.null(mode.prior.size)){
      if(mode.prior.size < n){mode.prior.size = n}
      beta <- 2*mode.prior.size/n - 1
      if(is.null(alpha)) alpha=mode.prior.size/(mode.prior.size-n)
      if(is.null(maxN)){maxN <- min(maxNmax,ceiling( n/(1-0.90^(1/beta)) ))}
      fn <- function(beta,x,n,mode.prior.size,effective.prior.df,alpha){
       priorm <- dfn(alpha,beta,x,n,effective.prior.df)
# print(c(beta,sum(x*priorm)/sum(priorm),
#         abs(mode.prior.size - (x+0.5)[which.max(priorm)]) ))
       abs(mode.prior.size - (x+0.5)[which.max(priorm)])
      }
      maxN = ceiling(3*mode.prior.size)
      x <- n:maxN
      a = optimize(f=fn,interval=c(1,maxbeta),x,n,mode.prior.size,
                   effective.prior.df,alpha,tol=0.01)
      beta <- a$minimum
      while( {
       p=dfn(alpha,beta,x,n,effective.prior.df);
#      print( c(sum(p[x>1000]), sum(p[x>4000])) );
       abs(p[length(p)]/max(p) - 0.01)>0.005}){
        maxN <- round(maxN*(c(0.9,1.1)[(p[length(p)]/max(p) > 0.01)+1]))
        x <- n:maxN
        a = optimize(f=fn,interval=c(1,maxbeta),x,n,mode.prior.size,
                     effective.prior.df,alpha,tol=0.01)
        beta <- a$minimum
#print(c(beta,maxN,pbeta(n/maxN,shape1=alpha,shape2=beta)))
      }
      maxN = min(maxNmax,maxN)
     }
     if(!is.null(quartiles.prior.size)){
      if(quartiles.prior.size[1] < n){quartiles.prior.size[1] = n}
      if(quartiles.prior.size[2] < quartiles.prior.size[1]){
        aaa <- quartiles.prior.size[2]
        quartiles.prior.size[2] = quartiles.prior.size[1]
        quartiles.prior.size[1] = aaa
      }
      fn <- function(p,x,n,quartiles.prior.size,effective.prior.df){
       priorm <- dfn(exp(p[1]),exp(p[2]),x,n,effective.prior.df)
# print(c(exp(p),sum(x*priorm)/sum(priorm),
#      sqrt((quartiles.prior.size[1] - (x+0.5)[match(TRUE,cumsum(priorm) >= 0.25)])^2+ 
#           (quartiles.prior.size[2] - (x+0.5)[match(TRUE,cumsum(priorm) >= 0.75)])^2)))
       sqrt((quartiles.prior.size[1] - (x+0.5)[match(TRUE,cumsum(priorm) >= 0.25)])^2+ 
            (quartiles.prior.size[2] - (x+0.5)[match(TRUE,cumsum(priorm) >= 0.75)])^2)
      }
      maxN = ceiling(10*quartiles.prior.size[2])
      x <- n:maxN
      a = optim(par=log(c(1,10)),fn=fn,
        x=x,n=n,quartiles.prior.size=quartiles.prior.size,effective.prior.df=effective.prior.df,
        control=list(abstol=10))
      alpha <- exp(a$par[1])
      beta  <- exp(a$par[2])
      while( {
       p=dfn(alpha,beta,x,n,effective.prior.df);
#      print( c(sum(p[x>1000]), sum(p[x>4000])) );
       abs(p[length(p)]/max(p) - 0.01)>0.005}){
        maxN <- round(maxN*(c(0.9,1.1)[(p[length(p)]/max(p) > 0.01)+1]))
        x <- n:maxN
        a = optim(par=a$par,fn=fn,
         x=x,n=n,quartiles.prior.size=quartiles.prior.size,effective.prior.df=effective.prior.df,
         control=list(abstol=10))
        alpha <- exp(a$par[1])
        beta  <- exp(a$par[2])
#print(c(alpha,beta,maxN,sum(p[x>quartiles.prior.size[1]]), sum(p[x>quartiles.prior.size[2]] )))
      }
      maxN = min(maxNmax,maxN)
     }
     }
     if(is.null(alpha)) alpha=1
     if(is.null(maxN)){maxN <- min(maxNmax,ceiling( n/(1-0.90^(1/beta)) ))}
     if(is.null(N)){N <- min(maxNmax,ceiling( n/(1-0.5^(1/beta)) ))}
     x <- n:maxN
     lpriorm=dfn(alpha,beta,x,n,effective.prior.df);
     if(log){
      lpriorm <- log(lpriorm)
     }
     lpriorm
     }
    )
#    End of switch to compute lpriorm
     x <- n:maxN
     if(log){
      priorm <- exp(lpriorm)
     }else{
      priorm <- lpriorm
     }
     mode.prior.size <- (x+0.5)[which.max(priorm)]
     mean.prior.size <- sum((x+0.5)*priorm)
     median.prior.size <- (x+0.5)[match(TRUE,cumsum(priorm) >= 0.5)]
     quartiles.prior.size[1] <- (x+0.5)[match(TRUE,cumsum(priorm) >= 0.25)]
     quartiles.prior.size[2] <- (x+0.5)[match(TRUE,cumsum(priorm) >= 0.75)]
    if(is.null(N)){N <- mean(x)}
    if(verbose){cat(paste("The maximum prior population size is",maxN,"\n"))}
    list(x=x,lprior=lpriorm,N=N,maxN=maxN,
         median.prior.size=median.prior.size,
         mean.prior.size=mean.prior.size,
         mode.prior.size=mode.prior.size,
         quartiles.prior.size=quartiles.prior.size,
	 mode.prior.sample.proportion=mode.prior.sample.proportion,
	 median.prior.sample.proportion=median.prior.sample.proportion,
	 alpha=alpha,beta=beta,
	 effective.prior.df=effective.prior.df,
         type=type)
}
