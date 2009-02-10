margposteriorsize<-function(s,N=trunc(length(s)*seq(1.1,4,length=10)+1),
          K=(1:max(s)),
          mu=5,rho=3,M=100000,
          parallel=1,seed=NULL, verbose=FALSE, 
          n=tabulate(s,nbin=K),
          lCval=NULL,
          nstart=NULL,
          return.all=TRUE){
  if(length(N)==0||!is.numeric(N)||trunc(N)!=N){
    stop("N needs to be a vector of integers.")
  }
  if(any(length(s)>N)){print("Error - observed unit outside range")}
  if(!is.null(seed))  set.seed(as.integer(seed))
  prob=rep(0,K)
  llik <- function(x,s,n){ llhoodf(N=exp(x),s=s,n=n) }
  if(is.null(lCval)){
   if(is.null(nstart)){
          nstart=n/(1:K)
          nstart=log(3*sum(n)*nstart/sum(nstart))
   }
   out <- optim(par=nstart, fn=llik, s=s, n=n,
    control=list(maxit=10000,fnscale=-1))
   Nmle <- exp(out$par)
   lCval <- out$value
   print(lCval)
  }else{
   Nmle <- nstart
  }
  if(parallel==1){
    Cret <- .C("bnw_mp",
              N=as.integer(N),
              lenN=as.integer(length(N)),
              K=as.integer(K),
              n=as.integer(length(s)),
              s=as.integer(s),
              nk=as.integer(n),
              lCval=as.double(lCval),
              prob=as.double(N),
              Nprior=as.integer(prob),
              Nmle=as.integer(round(Nmle)),
              mu=as.double(mu),
              rho=as.double(rho),
              M=as.integer(M))
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
      cat("no problem... PVM started by rpm...\n")
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
    outlist <- clusterCall(cl, margposteriorsize,
      s=s,N=N,K=K,mu=mu,rho=rho,n=n,lCval=lCval,nstart=Nmle,
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
    Cret$prob <- Cret$prob / np
    for(i in 1:np){
      out <- outlist[[i]]
      Cret$prob <- Cret$prob + out$prob / np
      if(out$lCval > Cret$lCval){
        Cret$lCval <- out$lCval
        Cret$Nmle <- out$Nmle
      }
    }
    if(verbose){
     cat("parallel samplesize=", parallel,"by", MCsamplesize.parallel,"\n")
    }
#   }
    stopCluster(cl)
    if(getClusterOption("type")=="PVM") .PVM.exit()
  }
  if(return.all){
    Cret
  }else{
    Cret$prob
  }
}
margpossize<-function(s,N=trunc(length(s)*seq(1.1,4,length=10)+1),
          K=(1:max(s)),
          mu=5,rho=3,M=100000,
          parallel=1,seed=NULL, verbose=FALSE, 
          n=tabulate(s,nbin=K),
          return.all=FALSE){
  if(length(N)==0||!is.numeric(N)||trunc(N)!=N){
    stop("N needs to be a vector of integers.")
  }
  logprob <- N
  for(i in seq(along=N)){
   logprob[i] <- margposN(N=N[i],K=K,s=s,mu=mu,rho=rho,M=M,
               parallel=parallel,seed=seed, verbose=verbose, 
               n=n, return.all=return.all)
   if(verbose) cat(paste("N=",N[i]," logprob=",logprob[i],"\n"))
  }
  list(N=N,logprob=logprob,s=s,K=K,mu=mu,rho=rho)
}
unposN<-function(mu,rho,N,K,s,n=tabulate(s,nbin=K),M=100000){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    if(length(s)>N){print("Error - observed unit outside range")}
    prob=rep(0,K)
    Cret <- .C("bnw_NC",
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
              unpos=as.double(0))
    Cret$unpos
    Cret
}
margposN<-function(N,K,s,mu=5,rho=3,M=100000,
          parallel=1,seed=NULL, verbose=FALSE, 
          n=tabulate(s,nbin=K),
          return.all=FALSE){
  #this function takes a vector of population sizes and a vector s of 
  #sequential sizes of sampled units and returns a log likelihood value
  #s values must all be positive integers
  if(length(s)>N){print("Error - observed unit outside range")}
  if(!is.null(seed))  set.seed(as.integer(seed))
  prob=rep(0,K)
  if(parallel==1){
    Cret <- .C("bnw_NC",
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
      cat("no problem... PVM started by rpm...\n")
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
      mu=mu,rho=rho,N=N,K=K,s=s,M=MCsamplesize.parallel)
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

llhoodf<-function(N,s,n=tabulate(s,nbin=length(N))){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    if(length(s)>sum(N)){print("Error - observed unit outside range")}
    Cret <- .C("bnw_llik",
              K=as.integer(length(N)),
              n=as.integer(length(s)),
              s=as.integer(s),
              snk=as.integer(n),
              Nk=as.double(N),
              llik=as.double(0))
    Cret$llik
}
