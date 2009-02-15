margposteriorsize<-function(N=trunc(length(s)*seq(1.1,4,length=10)+1),
          s,
	  K=max(s),
          prob=rep(1/K,K),
          M=100000,
          parallel=1,seed=NULL, verbose=TRUE, 
          n=tabulate(s,nbin=K),
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
    outlist <- clusterApply(cl, as.list(N), margposN,
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
    control=list(maxit=1000,fnscale=-1))
   Nmle <- toN(out$par, Ntot=N, nk=n)
#  lbound <- llhoodf(N=Nmle, s=s, verbose=FALSE)
#  qprob <- n/sum(n)
   qprob <- Nmle/sum(Nmle)
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
      s=s,N=N,K=K,prob=prob,M=MCsamplesize.parallel,n=n,qprob=qprob)
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
	  qprob=NULL,
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
    .C("bnw_llik",
              K=as.integer(length(N)),
              n=as.integer(length(s)),
              s=as.integer(s),
              nk=as.integer(n),
              Nk=as.double(N),
              llik=as.double(0))$llik
  }
  if(is.null(qprob)){
   out <- N*n/sum(n)
   out <- optim(par=tox(out,nk=n), fn=llik, 
     Ntot=N, n=n, s=s,
    control=list(maxit=10000,fnscale=-1))
   Nmle <- toN(out$par, Ntot=N, nk=n)
   lbound <- llhoodf(N=Nmle, s=s, verbose=FALSE)
   qprob <- Nmle/sum(Nmle)
  }
# print(qprob)
# print(lbound)
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
      s=s,N=N,K=K,M=MCsamplesize.parallel,n=n,qprob=qprob)
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
          parallel=1,seed=NULL, verbose=TRUE, 
          n=tabulate(s,nbin=K),
          return.all=TRUE){
  if(length(N)==0||!is.numeric(N)||trunc(N)!=N){
    stop("N needs to be a vector of integers.")
  }
  if(!is.null(seed))  set.seed(as.integer(seed))

  if(length(N)==1){
      Cret <- discretemleN(s=s,N=N,K=K,M=M,seed=seed,n=n,
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
    outlist <- clusterApply(cl, as.list(N), discretemleN,
      s=s,K=K,M=M,n=n,return.all=TRUE,parallel=1)
#
#   Process the results
#
    Cret <- outlist[[1]]
    Cret$N <- N
    Cret$llik <- Cret$mllik
    Cret$Nmle <- Cret$Nk
    Cret$Nk <- NULL
    for(i in (2 : length(outlist))){
     z <- outlist[[i]]
     Cret$llik <- c(Cret$llik,z$mllik)
     Cret$lbound <- c(Cret$lbound,z$lbound)
     if(z$mllik > Cret$mllik){
      Cret$Nmle <- z$Nk
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
