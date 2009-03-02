margposteriorsize<-function(N=trunc(length(s)*seq(1.1,4,length=10)+1),
          s,
	  K=max(s),
          prob=rep(1/K,K),
          M=100000,
          parallel=1, seed=NULL, verbose=TRUE, 
          n=tabulate(s,nbin=K),
	  maxit=50000, method="Nelder-Mead", temp=10, trace=0,
	  return.all=TRUE){
  if(length(prob)!=K){
    stop("The vector of probabilities need to have length K.")
  }
  if(length(N)==0||!is.numeric(N)||trunc(N)!=N){
    stop("N needs to be a vector of integers.")
  }
  if(!is.null(seed))  set.seed(as.integer(seed))

  if(length(N)==1){
      Cret <- margposN(s=s,N=N,K=K,prob=prob,M=M,seed=seed,n=n,
                       return.all=return.all,parallel=parallel,
		       maxit=maxit,method=method,trace=trace,temp=temp)
  }else{
    require(snow)
#
#   Start PVM if necessary
#
#   setDefaultClusterOptions(type="PVM")
    if(getClusterOption("type")=="PVM") {
     if(verbose)
     {
      cat("Engaging warp drive using PVM ...\n")
     }
     require(rpvm)
     PVM.running <- try(.PVM.config(), silent = TRUE)
     if(inherits(PVM.running, "try-error"))
     {
      hostfile <- paste(Sys.getenv("HOME"), "/.xpvm_hosts", sep = "")
      .PVM.start.pvmd(hostfile)
      cat("no problem... PVM started by size...\n")
     }
    }else{
     if(verbose)
     {
      cat("Engaging warp drive using MPI ...\n")
     }
    }
#
#   Start Cluster
#
    cl <- makeCluster(parallel)
    clusterSetupRNG(cl)
    clusterEvalQ(cl, library(size))
#
#   Run the jobs with rpvm or Rmpi
#
    flush.console()
    outlist <- clusterApplyLB(cl, as.list(N), margposN,
      s=s,K=K,prob=prob,M=M,n=n,return.all=FALSE,parallel=1)
#
#   Process the results
#
    Cret <- list(unpos=unlist(outlist), N=Nrange, s=s,prob=prob,
                 M=M,n=n)
    if(verbose){
     cat("parallel samplesize=", parallel,"by", M,"\n")
    }
    stopCluster(cl)
    if(getClusterOption("type")=="PVM") .PVM.exit()
  }
  if(return.all){
    Cret
  }else{
    Cret$unpos
  }
}
unposN<-function(prob,N,K,s,n=tabulate(s,nbin=K),M=100000){
  #this function takes a vector of population sizes and a vector s of 
  #sequential sizes of sampled units and returns a log likelihood value
  #s values must all be positive integers
  if(length(s)>N){print("Error: The population counts should not be exceeded by the sample counts.")}
  Cret <- .C("bnw_NC",
              N=as.integer(N),
              K=as.integer(K),
              n=as.integer(length(s)),
              s=as.integer(s),
              nk=as.integer(n),
              Nk=as.integer(prob),
              prob=as.double(prob),
              qprob=as.double(prob),
              M=as.integer(M),
              unpos=as.double(0))
    Cret$unpos
    Cret
}
margposN<-function(N, s,
          K=max(s),
          prob=rep(1/K,K),
          M=100000,
          parallel=1, seed=NULL, verbose=FALSE, 
          n=tabulate(s,nbin=K),
	  qprob=NULL,
	  maxit=50000, method="Nelder-Mead", temp=10, trace=0,
          return.all=FALSE){
  #this function takes a vector of population sizes and a vector s of 
  #sequential sizes of sampled units and returns a log likelihood value
  #s values must all be positive integers
  if(length(s)>N){print("Error: The population counts should not be exceeded by the sample counts.")}
  if(!is.null(seed))  set.seed(as.integer(seed))
  toN <- function(x,Ntot,nk){
    nk[nk!=0] <- (Ntot-sum(nk))*c(1,exp(x))/(1+sum(exp(x)))+nk[nk!=0]
    nk
  }
  tox <- function(N,nk){
    N <- N[nk!=0]-nk[nk!=0]
    ltox <- log(N[-1]/N[1])
    ltox[is.na(ltox) || is.infinite(ltox)] <- -100000
    ltox
  } 
  llik <- function(x,Ntot,n,s,prob,M=1000){
    qprob <- toN(x,Ntot,n)
    qprob <- qprob/sum(qprob)
    .C("bnw_NCbound",
       N=as.integer(Ntot),
       K=as.integer(length(qprob)),
       n=as.integer(length(s)),
       s=as.integer(s),
       nk=as.integer(n),
       Nk=as.integer(n),
       prob=as.double(prob),
       qprob=as.double(qprob),
       M=as.integer(M),
       unpos=as.double(0))$unpos
  }
  if(is.null(qprob)){
   out <- N*n/sum(n)
   out <- optim(par=tox(out,nk=n), fn=llik, 
     Ntot=N, n=n, s=s, prob=prob,
     method=method,
    control=list(maxit=maxit,fnscale=-1,temp=temp,trace=trace))
   Nkmle <- toN(out$par, Ntot=N, nk=n)
#  lbound <- llhoodf(N=Nkmle, s=s, verbose=FALSE)
#  qprob <- n/sum(n)
   qprob <- Nkmle/sum(Nkmle)
  }
  print(qprob)
  if(parallel==1){
    Cret <- .C("bnw_NC",
              N=as.integer(N),
              K=as.integer(K),
              n=as.integer(length(s)),
              s=as.integer(s),
              nk=as.integer(n),
              Nk=as.integer(n),
              prob=as.double(prob),
              qprob=as.double(qprob),
              M=as.integer(M),
              unpos=as.double(0))
  }else{
    require(snow)
#
#   Start PVM if necessary
#
#   setDefaultClusterOptions(type="PVM")
    if(getClusterOption("type")=="PVM") {
     if(verbose)
     {
      cat("Engaging warp drive using PVM ...\n")
     }
     require(rpvm)
     PVM.running <- try(.PVM.config(), silent = TRUE)
     if(inherits(PVM.running, "try-error"))
     {
      hostfile <- paste(Sys.getenv("HOME"), "/.xpvm_hosts", sep = "")
      .PVM.start.pvmd(hostfile)
      cat("no problem... PVM started by size...\n")
     }
    }else{
     if(verbose)
     {
      cat("Engaging warp drive using MPI ...\n")
     }
    }
#
#   Start Cluster
#
    cl <- makeCluster(parallel)
    clusterSetupRNG(cl)
    clusterEvalQ(cl, library(size))
#
#   Run the jobs with rpvm or Rmpi
#
    flush.console()
    MCsamplesize.parallel=round(M/parallel)
    outlist <- clusterCall(cl, margposN,
      s=s,N=N,K=K,prob=prob,M=MCsamplesize.parallel,n=n,qprob=qprob,
      maxit=maxit,method=method,trace=trace,temp=temp)
#
#   Process the results
#
#   Cret <- 0
#   for(i in (1 : parallel)){
#    z <- outlist[[i]]
#    Cret <- Cret+z/parallel
#   }
    Cret <- list(unpos=mean(unlist(outlist)))
    if(verbose){
     cat("parallel samplesize=", parallel,"by", MCsamplesize.parallel,"\n")
    }
    stopCluster(cl)
    if(getClusterOption("type")=="PVM") .PVM.exit()
  }
  Cret$Nmle <- sum(Cret$Nk)
  if(return.all){
    Cret
  }else{
    Cret$unpos
  }
}

llhoodf<-function(N,s,n=tabulate(s,nbin=length(N)),verbose=TRUE){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    if(any(n>N) || length(s)>sum(N)){
      if(verbose) print("Error: The population counts should not be exceeded by the sample counts.")
      return(-1000000)
    }
    Cret <- .C("bnw_llik",
              K=as.integer(length(N)),
              n=as.integer(length(s)),
              s=as.integer(s),
              nk=as.integer(n),
              Nk=as.double(N),
              llik=as.double(0))
    Cret$llik
}
discretemleN<-function(N,s,
          K=max(s),
          M=100000,
          parallel=1, seed=NULL, verbose=FALSE, 
          n=tabulate(s,nbin=K),
          Nkmle.start=NULL,
	  qprob=NULL,
	  maxit=50000, maxit.spps=10, method="Nelder-Mead", temp=10, trace=0,
          return.all=TRUE){
  if(length(s)>N){print("Error: The population counts should not be exceeded by the sample counts.")}
  if(!is.null(seed))  set.seed(as.integer(seed))
  toN <- function(x,Ntot,nk){
    nk[nk!=0] <- (Ntot-sum(nk))*c(1,exp(x))/(1+sum(exp(x)))+nk[nk!=0]
    nk
  }
  tox <- function(N,nk){
    N <- N[nk!=0]-nk[nk!=0]
    ltox <- log(N[-1]/N[1])
    ltox[is.na(ltox) || is.infinite(ltox)] <- -100000
    ltox
  } 
  llik <- function(x,Ntot,n,s){
    N <- toN(x,Ntot,n)
    temp <- .C("bnw_llik",
        K=as.integer(length(N)),
        n=as.integer(length(s)),
        s=as.integer(s),
        nk=as.integer(n),
        Nk=as.double(N),
        llik=as.double(0))$llik
    if(is.na(temp) || is.infinite(temp)) {temp <- -100000}
    temp
  }
  if(is.null(qprob)){
   if(is.null(Nkmle.start)){
     Nkmle.start <- sppsestN(N=N,s=s,maxit=maxit.spps)
     Nkmle.start <- abs(Nkmle.start-n)
     Nkmle.start <- n+(N-sum(n))*Nkmle.start/sum(Nkmle.start)
#    print(cbind(Nkmle.start,n))
#    if(any(Nkmle.start<n)){
#      Nkmle.start <- akestN(N=N,s=s)
#      Nkmle.start <- abs(Nkmle.start-n)
#      Nkmle.start <- n+(N-sum(n))*Nkmle.start/sum(Nkmle.start)
#    }
     if(any(Nkmle.start<n)){
       print(paste("Nkmle.start=",Nkmle.start))
       print(paste("nk=",n))
       print(cbind(Nkmle.start,n))
       Nkmle.start[Nkmle.start<n] <- n[Nkmle.start<n]+0.1
       if(Nkmle.start[1]==0){Nkmle.start[1]=0.1}
       Nkmle.start <- N*Nkmle.start/sum(Nkmle.start)
     }
   }
   out <- optim(par=tox(Nkmle.start,nk=n), fn=llik, 
     Ntot=N, n=n, s=s,
     method=method,
    control=list(maxit=maxit,fnscale=-1,trace=trace,temp=temp))
   Nkmle <- toN(out$par, Ntot=N, nk=n)
   lbound <- llhoodf(N=Nkmle, s=s, verbose=FALSE)
#  if(any(Nkmle<n)){
#      print(paste("Nkmle=",Nkmle))
#      print(paste("nk=",n))
#      Nkmle[Nkmle<n] <- n[Nkmle<n]
#  }
   Nkmle <- abs(Nkmle-n)
#  Nkmle <- n+(N-sum(n))*Nkmle/sum(Nkmle)
   qprob <- Nkmle/sum(Nkmle)
  }else{
   lbound <- -1000000
  }
  if(parallel==1){
    Cret <- .C("bnw_stocdiscrete",
              N=as.integer(N),
              K=as.integer(K),
              n=as.integer(length(s)),
              s=as.integer(s),
              nk=as.integer(n),
              Nk=as.integer(n),
              qprob=as.double(qprob),
              M=as.integer(M),
              mllik=as.double(0))
    Cret$llik <- Cret$mllik
    Cret$Nkmle <- Cret$Nk
    Cret$lbound <- lbound
  }else{
    require(snow)
#
#   Start PVM if necessary
#
#   setDefaultClusterOptions(type="PVM")
    if(getClusterOption("type")=="PVM") {
     if(verbose)
     {
      cat("Engaging warp drive using PVM ...\n")
     }
     require(rpvm)
     PVM.running <- try(.PVM.config(), silent = TRUE)
     if(inherits(PVM.running, "try-error"))
     {
      hostfile <- paste(Sys.getenv("HOME"), "/.xpvm_hosts", sep = "")
      .PVM.start.pvmd(hostfile)
      cat("no problem... PVM started by size...\n")
     }
    }else{
     if(verbose)
     {
      cat("Engaging warp drive using MPI ...\n")
     }
    }
#
#   Start Cluster
#
    cl <- makeCluster(parallel)
    clusterSetupRNG(cl)
    clusterEvalQ(cl, library(size))
#
#   Run the jobs with rpvm or Rmpi
#
    flush.console()
    MCsamplesize.parallel=round(M/parallel)
    outlist <- clusterCall(cl, discretemleN,
      s=s,N=N,K=K,M=MCsamplesize.parallel,n=n,Nkmle.start=Nkmle.start,
      qprob=qprob, maxit=maxit,method=method,trace=trace,temp=temp)
#
#   Process the results
#
    Cret <- outlist[[1]]
    for(i in (2 : length(outlist))){
     z <- outlist[[i]]
     if(z$mllik < Cret$mllik){
      Cret$Nk <- z$Nk
      Cret$mllik <- z$mllik
     }
    }
    if(verbose){
     cat("parallel samplesize=", parallel,"by", MCsamplesize.parallel,"\n")
    }
    stopCluster(cl)
    if(getClusterOption("type")=="PVM") .PVM.exit()
  }
  Cret$lbound <- lbound
  Cret$Nmle <- sum(Cret$Nk)
  if(return.all){
    Cret
  }else{
    Cret$mllik
  }
}
discretemleNimpute<-function(N,s,
          nsim,
          K=max(s),
          M=100000,
          parallel=1, seed=NULL, verbose=FALSE, 
	  n=NULL,
	  qprob=NULL,
	  maxit=50000, method="Nelder-Mead", temp=10, trace=0,
          return.all=TRUE){
  if(length(s)/nsim>N){print("Error: The population counts should not be exceeded by the sample counts.")}
  if(!is.null(seed))  set.seed(as.integer(seed))
  ns <- length(s)/nsim
  if(is.null(n)){
   n <- as.vector(apply(matrix(s,nrow=nsim,byrow=T),1,tabulate,nbin=K))
  }
  toN <- function(x,Ntot,nk){
    nk[nk!=0] <- (Ntot-sum(nk))*c(1,exp(x))/(1+sum(exp(x)))+nk[nk!=0]
    nk
  }
  tox <- function(N,nk){
    N <- N[nk!=0]-nk[nk!=0]
    ltox <- log(N[-1]/N[1])
    ltox[is.na(ltox) || is.infinite(ltox)] <- -100000
    ltox
  } 
  llik <- function(x,Ntot,n,s){
    N <- toN(x,Ntot,n)
    .C("bnw_llik",
              K=as.integer(length(N)),
              n=as.integer(length(s)),
              s=as.integer(s),
              nk=as.integer(n),
              Nk=as.double(N),
              llik=as.double(0))$llik
  }
  if(is.null(qprob)){
   cprob <- NULL
   lbound <- NULL
   for(i in 1:nsim){
    ni <- n[(i-1)*K+(1:K)]
    si <- s[(i-1)*ns+(1:ns)]
    out <- N*ni/sum(ni)
    out <- optim(par=tox(out,nk=ni), fn=llik, 
     Ntot=N, n=ni, s=si,
     method=method,
     control=list(maxit=maxit,fnscale=-1,trace=trace,temp=temp))
    Nkmle <- toN(out$par, Ntot=N, nk=ni)
    lbound <- c(lbound,llhoodf(N=Nkmle, s=si, verbose=FALSE))
    qprob <- c(qprob,Nkmle/sum(Nkmle))
   }
  }
  if(parallel==1){
    Cret <- .C("bnw_stocdiscreteimpute",
              N=as.integer(N),
              K=as.integer(K),
              n=as.integer(ns),
              s=as.integer(s),
              nk=as.integer(n),
              Nk=as.integer(n),
              qprob=as.double(qprob),
              nsim=as.integer(nsim),
              M=as.integer(M),
              mllik=as.double(0))
  }else{
    require(snow)
#
#   Start PVM if necessary
#
#   setDefaultClusterOptions(type="PVM")
    if(getClusterOption("type")=="PVM") {
     if(verbose)
     {
      cat("Engaging warp drive using PVM ...\n")
     }
     require(rpvm)
     PVM.running <- try(.PVM.config(), silent = TRUE)
     if(inherits(PVM.running, "try-error"))
     {
      hostfile <- paste(Sys.getenv("HOME"), "/.xpvm_hosts", sep = "")
      .PVM.start.pvmd(hostfile)
      cat("no problem... PVM started by size...\n")
     }
    }else{
     if(verbose)
     {
      cat("Engaging warp drive using MPI ...\n")
     }
    }
#
#   Start Cluster
#
    cl <- makeCluster(parallel)
    clusterSetupRNG(cl)
    clusterEvalQ(cl, library(size))
#
#   Run the jobs with rpvm or Rmpi
#
    flush.console()
    MCsamplesize.parallel=round(M/parallel)
    outlist <- clusterCall(cl, discretemleNimpute,
      s=s,N=N,nsim=nsim,K=K,M=MCsamplesize.parallel,qprob=qprob,
      maxit=maxit,method=method,trace=trace,temp=temp)
#
#   Process the results
#
    Cret <- outlist[[1]]
    for(i in (2 : length(outlist))){
     z <- outlist[[i]]
     if(z$mllik < Cret$mllik){
      Cret$Nk <- z$Nk
      Cret$mllik <- z$mllik
     }
    }
    if(verbose){
     cat("parallel samplesize=", parallel,"by", MCsamplesize.parallel,"\n")
    }
    stopCluster(cl)
    if(getClusterOption("type")=="PVM") .PVM.exit()
  }
  Cret$Nmle <- sum(Cret$Nkmle)
  if(return.all){
    Cret
  }else{
    Cret$mllik
  }
}
discretemle<-function(N=trunc(length(s)*seq(1.1,4,length=10)+1),
          s,
	  K=max(s),
          M=100000,
          parallel=1, seed=NULL, verbose=FALSE, 
          n=tabulate(s,nbin=K),
          Nkmle.start=NULL,
	  maxit=50000, method="Nelder-Mead", temp=10, trace=0,
          return.all=TRUE){
  if(length(N)==0||!is.numeric(N)||trunc(N)!=N){
    stop("N needs to be a vector of integers.")
  }
  if(!is.null(seed))  set.seed(as.integer(seed))

  if(length(N)==1){
      Cret <- discretemleN(s=s,N=N,K=K,M=M,seed=seed,n=n,
                Nkmle.start=Nkmle.start,maxit=maxit,method=method,trace=trace,temp=temp,
		return.all=return.all,parallel=parallel)
  }else{
    require(snow)
#
#   Start PVM if necessary
#
#   setDefaultClusterOptions(type="PVM")
    if(getClusterOption("type")=="PVM") {
     if(verbose)
     {
      cat("Engaging warp drive using PVM ...\n")
     }
     require(rpvm)
     PVM.running <- try(.PVM.config(), silent = TRUE)
     if(inherits(PVM.running, "try-error"))
     {
      hostfile <- paste(Sys.getenv("HOME"), "/.xpvm_hosts", sep = "")
      .PVM.start.pvmd(hostfile)
      cat("no problem... PVM started by size...\n")
     }
    }else{
     if(verbose)
     {
      cat("Engaging warp drive using MPI ...\n")
     }
    }
#
#   Start Cluster
#
    cl <- makeCluster(parallel)
    clusterSetupRNG(cl)
    clusterEvalQ(cl, library(size))
#
#   Run the jobs with rpvm or Rmpi
#
    flush.console()
    outlist <- clusterApplyLB(cl, as.list(N), discretemleN,
      s=s,K=K,M=M,n=n,Nkmle.start=Nkmle.start,return.all=TRUE,parallel=1,maxit=maxit,method=method,
      trace=trace,temp=temp)
#
#   Process the results
#
    Cret <- outlist[[1]]
    Cret$N <- N
    Cret$llik <- Cret$mllik
    Cret$Nkmle <- Cret$Nk
    Cret$Nk <- NULL
    for(i in (2 : length(outlist))){
     z <- outlist[[i]]
     Cret$llik <- c(Cret$llik,z$mllik)
     Cret$lbound <- c(Cret$lbound,z$lbound)
     if(z$mllik > Cret$mllik){
      Cret$Nkmle <- z$Nk
      Cret$mllik <- z$mllik
     }
    }
    stopCluster(cl)
    if(getClusterOption("type")=="PVM") .PVM.exit()
  }
  Cret$Nmle <- sum(Cret$Nmle)
  if(return.all){
    Cret
  }else{
    Cret$mllik
  }
}
discretemlesearch<-function(s,
          interval=trunc(length(s)*seq(1.1,4,length=10)+1),
	  K=max(s),
          M=10000,
          parallel=1, seed=NULL, verbose=FALSE, 
          n=tabulate(s,nbin=K),
          Nkmle.start=NULL,
	  maxit=10000, tol=0.1, maxit.spps=3,
          return.all=TRUE){
  if(!is.null(seed))  set.seed(as.integer(seed))

  llik <- function(x,s,K,M,seed,n,parallel,maxit,maxit.spps){
    discretemleN(N=round(x),s=s,K=K,M=M,seed=seed,n=n,
	     maxit=maxit,return.all=FALSE,parallel=parallel,maxit.spps=maxit.spps)
  }
#  if(is.null(Nkmle.start)){
#    Nkmle.start <- sppsestN(N=N,s=s)
#    Nkmle.start <- abs(Nkmle.start-n)
#    Nkmle.start <- n+(N-sum(n))*Nkmle.start/sum(Nkmle.start)
#    if(any(Nkmle.start<n)){
#      print(paste("Nkmle.start=",Nkmle.start))
#      print(paste("nk=",n))
#      print(cbind(Nkmle.start,n))
#      Nkmle.start[Nkmle.start<n] <- n[Nkmle.start<n]+0.1
#      if(Nkmle.start[1]==0){Nkmle.start[1]=0.1}
#      Nkmle.start <- N*Nkmle.start/sum(Nkmle.start)
#    }
#  }
   out <- optimize(f=llik, 
      s=s,K=K,M=M,seed=seed,n=n,parallel=1,maxit=maxit,maxit.spps=maxit.spps,
      interval=interval, 
      maximum = TRUE,tol=tol)
  #
  if(return.all){
    out
  }else{
    out$maximum
  }
}
discretemleimpute<-function(N=trunc(length(s)*seq(1.1,4,length=10)+1),
          s,
	  K=max(s),
          M=100000,
          parallel=1, seed=NULL, verbose=FALSE, 
	  nsim=10, truncate.size=K,
	  maxit=50000, method="Nelder-Mead", temp=10, trace=0,
          return.all=TRUE){
  if(length(N)==0||!is.numeric(N)||trunc(N)!=N){
    stop("N needs to be a vector of integers.")
  }
  if(!is.null(seed))  set.seed(as.integer(seed))

  sizeN <- 1:K
  s.raw <- s
  s <- NULL
  n <- NULL
  for(i in 1:nsim){
   si <- round(impute.size(s.raw))
   si[si > truncate.size] <- truncate.size
   s <- c(s,si)
   n <- c(n,tabulate(si,nbin=max(sizeN))[sizeN])
  }
  if(length(N)==1){
      Cret <- discretemleNimpute(N=N,s=s,nsim=nsim,n=n,K=K,M=M,seed=seed,
                maxit=maxit,method=method,trace=trace,temp=temp,
		return.all=return.all,parallel=parallel)
  }else{
    require(snow)
#
#   Start PVM if necessary
#
#   setDefaultClusterOptions(type="PVM")
    if(getClusterOption("type")=="PVM") {
     if(verbose)
     {
      cat("Engaging warp drive using PVM ...\n")
     }
     require(rpvm)
     PVM.running <- try(.PVM.config(), silent = TRUE)
     if(inherits(PVM.running, "try-error"))
     {
      hostfile <- paste(Sys.getenv("HOME"), "/.xpvm_hosts", sep = "")
      .PVM.start.pvmd(hostfile)
      cat("no problem... PVM started by size...\n")
     }
    }else{
     if(verbose)
     {
      cat("Engaging warp drive using MPI ...\n")
     }
    }
#
#   Start Cluster
#
    cl <- makeCluster(parallel)
    clusterSetupRNG(cl)
    clusterEvalQ(cl, library(size))
#
#   Run the jobs with rpvm or Rmpi
#
    flush.console()
    outlist <- clusterApplyLB(cl, as.list(N), discretemleNimpute,
      s=s,nsim=nsim,K=K,M=M,n=n,return.all=TRUE,parallel=1,
      maxit=maxit,method=method,trace=trace,temp=temp)
#
#   Process the results
#
    Cret <- outlist[[1]]
    Cret$N <- N
    Cret$llik <- Cret$mllik
    Cret$Nkmle <- Cret$Nk
    Cret$Nk <- NULL
    for(i in (2 : length(outlist))){
     z <- outlist[[i]]
     Cret$llik <- c(Cret$llik,z$mllik)
     Cret$lbound <- c(Cret$lbound,z$lbound)
     if(z$mllik > Cret$mllik){
      Cret$Nkmle <- z$Nk
      Cret$mllik <- z$mllik
     }
    }
    stopCluster(cl)
    if(getClusterOption("type")=="PVM") .PVM.exit()
  }
  if(return.all){
    Cret
  }else{
    Cret$mllik
  }
}
discretemleimputeindividual<-function(N=trunc(length(s)*seq(1.1,4,length=10)+1),
          s,
	  K=max(s),
          M=100000,
          parallel=1, seed=NULL, verbose=FALSE, 
	  maxit=50000, method="Nelder-Mead", temp=10, trace=0,
	  nsim=10, truncate.size=K,
          return.all=TRUE){
  if(length(N)==0||!is.numeric(N)||trunc(N)!=N){
    stop("N needs to be a vector of integers.")
  }
  if(!is.null(seed))  set.seed(as.integer(seed))

  sizeN <- 1:K
  s.raw <- s
  outraw <- matrix(0,nrow=nsim,ncol=length(N))
  colnames(outraw) <- N
  rownames(outraw) <- paste("sim",1:nsim,sep="")
  outmle <- matrix(0,nrow=nsim,ncol=K)
  colnames(outmle) <- 1:K
  rownames(outmle) <- paste("sim",1:nsim,sep="")
  for(i in 1:nsim){
   s <- round(impute.size(s.raw))
   s[s > truncate.size] <- truncate.size
   stab <- tabulate(s,nbin=max(sizeN))[sizeN]
   out <- discretemle(N=N,s=s,K=K,M=M,seed=seed,n=stab,
             maxit=maxit,method=method,trace=trace,temp=temp,
	     verbose=verbose,
             return.all=return.all,parallel=parallel)
   outraw[i,] <- out$llik
   outmle[i,] <- out$Nkmle
  }
  list(llik=outraw, Nkmle=outmle, Nmle=sum(outmle))
}
