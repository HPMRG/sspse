margposteriorsize.origandgood<-function(s,N=trunc(length(s)*seq(1.1,4,length=10)+1),
          K=max(s),
          prob=rep(1/K,K),
          M=100000,
          parallel=1,seed=NULL, verbose=TRUE, 
          n=tabulate(s,nbins=K),
          lbound=NULL,
          nstart=NULL,
          return.all=TRUE){
  if(length(prob)!=K){
    stop("The vector of probabilities need to have length K.")
  }
  if(length(N)==0||!is.numeric(N)||trunc(N)!=N){
    stop("N needs to be a vector of integers.")
  }
  if(any(n>nstart)){print("Error: The population counts should not be exceeded by the sample counts.")}
  if(!is.null(seed))  set.seed(as.integer(seed))
  if(is.null(lbound)){
   toN <- function(x,Ntot,nk){
     nk[nk!=0] <- Ntot*c(1,exp(x))/(1+sum(exp(x)))
     nk
   }
   tox <- function(N,nk){
     N <- N[nk!=0]
     ltox <- log(N[-1]/N[1])
     ltox[is.na(ltox) || is.infinite(ltox)] <- -100000
     ltox
   } 
   llik <- function(x,Ntot,s,nk){
     llhoodf(N=round(toN(x,Ntot,nk)),s=s,verbose=FALSE)
   }
   if(is.null(nstart)){
          nstart=n/(1:K)
          nstart=3*length(s)*nstart/sum(nstart)
   }
   Nest <- nstart
   for(i in seq(along=N)){
     out <- N[i]*Nest/sum(Nest)
     out <- optim(par=tox(out,nk=n), fn=llik, Ntot=N[i], s=s, nk=n,
      control=list(maxit=10000,fnscale=-1))
#    prob <- toN(out$par, Ntot=N[i], nk=n)
#    Nest <- round(prob)
#    prob <- prob/sum(prob)
     Nest <- round(toN(out$par, Ntot=N[i],nk=n))
     out$value <- llhoodf(N=Nest, s=s, verbose=FALSE)
     if(out$value > lbound || i ==1){
       lbound <- out$value
       Nmle <- Nest
     }
   }
  }else{
   Nmle <- nstart
  }
  sprob <- n/sum(n)
  sprob <- Nmle/sum(Nmle)
  if(verbose) cat(paste("Rejection sampling bound set to",lbound,"\n"))
  if(parallel==1){
    Cret <- .C("bnw_mp",
              N=as.integer(N),
              lenN=as.integer(length(N)),
              K=as.integer(K),
              n=as.integer(length(s)),
              s=as.integer(s),
              nk=as.integer(n),
              lbound=as.double(lbound),
              dprob=as.double(sprob),
              prob=as.double(N),
              NtotMLE=as.double(N),
              Nprior=as.integer(prob),
              Nmle=as.integer(round(Nmle)),
              M=as.integer(M), PACKAGE="size")
	      Cret$prob <- Cret$prob*prob/sprob
  }else{
    cl <- beginsnow(parallel)
    MCsamplesize.parallel=round(M/parallel)
    outlist <- clusterCall(cl, margposteriorsize,
      s=s,N=N,K=K,n=n,prob=prob,lbound=lbound,nstart=Nmle,
      M=MCsamplesize.parallel)
#   if(length(N)==0||!is.numeric(N)||trunc(N)!=N){
#
#   Process the results
#
#   Cret <- 0
#   for(i in (1 : parallel)){
#    z <- outlist[[i]]
#    Cret <- Cret+z/parallel
#   }
    np <- length(outlist)
    Cret <- outlist[[1]]
    Cret$M <- M
    Cret$prob <- Cret$prob / np
    for(i in 1:np){
      out <- outlist[[i]]
      Cret$prob <- Cret$prob + out$prob / np
      if(out$lbound > Cret$lbound){
        Cret$lbound <- out$lbound
        Cret$Nmle <- out$Nmle
      }
    }
    Cret$prob <- Cret$prob*prob/sprob
    if(verbose){
     cat("parallel samplesize=", parallel,"by", MCsamplesize.parallel,"\n")
    }
    endsnow(cl)
  }
  if(return.all){
    Cret
  }else{
    Cret$prob
  }
}
margposteriorsizewar<-function(s,N=trunc(length(s)*seq(1.1,4,length=10)+1),
          K=max(s),
          mu=5,rho=3,M=100000,
          parallel=1,seed=NULL, verbose=TRUE, 
          n=tabulate(s,nbins=K),
          lbound=NULL,
          nstart=NULL,
          return.all=TRUE){
  if(length(N)==0||!is.numeric(N)||trunc(N)!=N){
    stop("N needs to be a vector of integers.")
  }
  if(any(n>nstart)){print("Error: The population counts should not be exceeded by the sample counts.")}
  if(!is.null(seed))  set.seed(as.integer(seed))
  if(is.null(lbound)){
   toN <- function(x,Ntot,nk){
     nk[nk!=0] <- Ntot*c(1,exp(x))/(1+sum(exp(x)))
     nk
   }
   tox <- function(N,nk){
     N <- N[nk!=0]
     ltox <- log(N[-1]/N[1])
     ltox[is.na(ltox) || is.infinite(ltox)] <- -100000
     ltox
   } 
   llik <- function(x,Ntot,s,nk){
     llhoodf(N=round(toN(x,Ntot,nk)),s=s,verbose=FALSE)
   }
   if(is.null(nstart)){
          nstart=n/(1:K)
          nstart=3*length(s)*nstart/sum(nstart)
   }
   for(i in seq(along=N)){
     out <- N[i]*nstart/sum(nstart)
     out <- optim(par=tox(out,nk=n), fn=llik, Ntot=N[i], s=s, nk=n,
      control=list(maxit=10000,fnscale=-1))
     if(out$value > lbound || i ==1){
       lbound <- out$value
     }
     if(i==1){
       Nmle <- round(toN(out$par, Ntot=N[i],nk=n))
     }
   }
  }else{
   Nmle <- nstart
  }
  print(Nmle)
  if(verbose) cat(paste("Rejection sampling bound set to",lbound,"\n"))
  if(parallel==1){
    Cret <- .C("bnw_mpwar",
              N=as.integer(N),
              lenN=as.integer(length(N)),
              K=as.integer(K),
              n=as.integer(length(s)),
              s=as.integer(s),
              nk=as.integer(n),
              lbound=as.double(lbound),
              prob=as.double(N),
              NtotMLE=as.double(N),
              Nprior=as.integer(prob),
              Nmle=as.integer(round(Nmle)),
              mu=as.double(mu),
              rho=as.double(rho),
              M=as.integer(M), PACKAGE="size")
  }else{
    cl <- beginsnow(parallel)
    MCsamplesize.parallel=round(M/parallel)
    outlist <- clusterCall(cl, margposteriorsizewar,
      s=s,N=N,K=K,mu=mu,rho=rho,n=n,lbound=lbound,nstart=Nmle,
      M=MCsamplesize.parallel)
#   if(length(N)==0||!is.numeric(N)||trunc(N)!=N){
#
#   Process the results
#
#   Cret <- 0
#   for(i in (1 : parallel)){
#    z <- outlist[[i]]
#    Cret <- Cret+z/parallel
#   }
    np <- length(outlist)
    Cret <- outlist[[1]]
    Cret$M <- M
    Cret$prob <- Cret$prob / np
    for(i in 1:np){
      out <- outlist[[i]]
      Cret$prob <- Cret$prob + out$prob / np
      if(out$lbound > Cret$lbound){
        Cret$lbound <- out$lbound
        Cret$Nmle <- out$Nmle
      }
    }
    if(verbose){
     cat("parallel samplesize=", parallel,"by", MCsamplesize.parallel,"\n")
    }
    endsnow(cl)
  }
  if(return.all){
    Cret
  }else{
    Cret$prob
  }
}
margpossize<-function(s,N=trunc(length(s)*seq(1.1,4,length=10)+1),
          K=max(s),
          M=100000,
          parallel=1,seed=NULL, verbose=FALSE, 
          n=tabulate(s,nbins=K),
          return.all=FALSE){
  if(length(N)==0||!is.numeric(N)||trunc(N)!=N){
    stop("N needs to be a vector of integers.")
  }
  logprob <- N
  for(i in seq(along=N)){
   logprob[i] <- margposN(N=N[i],K=K,s=s,M=M, #mu=mu,rho=rho,
               parallel=parallel,seed=seed, verbose=verbose, 
               n=n, return.all=return.all)
   if(verbose) cat(paste("N=",N[i]," logprob=",logprob[i],"\n"))
  }
  list(N=N,logprob=logprob,s=s,K=K,mu=mu,rho=rho)
}
unposNwar<-function(mu,rho,N,K,s,n=tabulate(s,nbins=K),M=100000){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    if(any(length(s)>N)){print("Error: The population counts should not be exceeded by the sample counts.")}
    prob=rep(0,K)
    Cret <- .C("bnw_NCwar",
              N=as.integer(N),
              K=as.integer(K),
              n=as.integer(length(s)),
              s=as.integer(s),
              nk=as.integer(n),
              Nk=as.integer(prob),
              prob=as.double(prob),
              mu=as.double(mu),
              rho=as.double(rho),
              M=as.integer(M),
              unpos=as.double(0), PACKAGE="size")
    Cret$unpos
    Cret
}
