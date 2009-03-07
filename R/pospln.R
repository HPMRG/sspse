poswest<-function(s,maxN=4*length(s),
                  K=2*max(s), nk=tabulate(s,nbin=K), n=length(s),
		  N=maxN/2,
                  prior.mean.degree=7, prior.sd.degree=2.2,
                  prior.mean.df=4, prior.sd.df=10,
                  muproposal=0.01, 
                  sigmaproposal=0.1, 
                  samplesize=10,burnin=0,interval=1,burnintheta=500,
                  seed=NULL,
                  verbose=TRUE){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    sigma0 <- sqrt(log(1+prior.sd.degree*prior.sd.degree/(prior.mean.degree*prior.mean.degree)))
    mu0 <- log(prior.mean.degree)-0.5*sigma0*sigma0
    if(!is.null(seed))  set.seed(as.integer(seed))
    Cret <- .C("gspps",
              pop=as.integer(c(s,rep(0,maxN-n))),
              nk=as.integer(nk),
              K=as.integer(K),
              n=as.integer(n),
              samplesize=as.integer(samplesize),
              burnin=as.integer(burnin),
              interval=as.integer(interval),
              mu0=as.double(mu0), prior.mean.df=as.double(prior.mean.df),
              sigma0=as.double(sigma0), prior.sd.df=as.double(prior.sd.df),
              muproposal=as.double(muproposal),
              sigmaproposal=as.double(sigmaproposal),
              N=as.integer(N),
              maxN=as.integer(maxN),
              sample=double(samplesize*4),
              burnintheta=as.integer(burnintheta),
              fVerbose=as.integer(verbose))
    Cret$sample<-matrix(Cret$sample,nrow=samplesize,ncol=4,byrow=TRUE)
    colnames(Cret$sample) <- c("N","mu","sigma","isolates")
    Cret$sample[,"mu"] <- exp(Cret$sample[,"mu"]+0.5*Cret$sample[,"sigma"]*Cret$sample[,"sigma"])
    Cret$sample[,"sigma"] <- Cret$sample[,"mu"]*sqrt(exp(Cret$sample[,"sigma"]*Cret$sample[,"sigma"])-1)
    Cret$nk<-Cret$nk/sum(Cret$nk)
    endrun <- burnin+interval*(samplesize-1)
    attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    mapfn <- function(x){
      posdensN <- density(x)
      posdensN$x[which.max(posdensN$y)]
    }
    Cret$MAP <- apply(Cret$sample,2,mapfn)
    Cret
}
poswestN<-function(s,N,
                 K=max(s), nk=tabulate(s,nbin=K), n=length(s),
                 prior.mean.degree=7, prior.sd.degree=2.2,
                 prior.mean.df=4,prior.sd.df=10,
                 muproposal=0.01, 
                 sigmaproposal=0.1, 
                 samplesize=10,burnin=0,interval=1,burnintheta=500,
                 seed=NULL,
                 verbose=TRUE){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    sigma0 <- sqrt(log(1+prior.sd.degree*prior.sd.degree/(prior.mean.degree*prior.mean.degree)))
    mu0 <- log(prior.mean.degree)-0.5*sigma0*sigma0
    if(!is.null(seed))  set.seed(as.integer(seed))
    Cret <- .C("gsppsN",
              pop=as.integer(c(s,rep(0,N-n))),
              nk=as.integer(nk),
              K=as.integer(K),
              n=as.integer(n),
              samplesize=as.integer(samplesize),
              burnin=as.integer(burnin),
              interval=as.integer(interval),
              mu0=as.double(mu0), prior.mean.df=as.double(prior.mean.df),
              sigma0=as.double(sigma0), prior.sd.df=as.double(prior.sd.df),
              muproposal=as.double(muproposal),
              sigmaproposal=as.double(sigmaproposal),
              N=as.integer(N),
              sample=integer(samplesize*N),
              burnintheta=as.integer(burnintheta),
              fVerbose=as.integer(verbose))
    Cret$sample<-matrix(Cret$sample,nrow=samplesize,ncol=N,byrow=TRUE)
    colnames(Cret$sample) <- paste(1:N)
    endrun <- burnin+interval*(samplesize-1)
    attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    Cret
}
pospln<-function(pop, K=max(pop),
                 prior.mean.degree=7, prior.sd.degree=2.2,
                 prior.mean.df=4,prior.sd.df=10,
                 muproposal=0.01, 
                 sigmaproposal=0.1, 
                 samplesize=1000,burnin=0,interval=1,
                 seed=NULL,
                 verbose=TRUE){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    sigma0 <- sqrt(log(1+prior.sd.degree*prior.sd.degree/(prior.mean.degree*prior.mean.degree)))
    mu0 <- log(prior.mean.degree)-0.5*sigma0*sigma0
    if(!is.null(seed))  set.seed(as.integer(seed))
    nk <- tabulate(pop,nbin=K)
    musample <- rep(0,samplesize)
    sigmasample <- rep(0,samplesize)
    Cret <- .C("MHpln",
              nk=as.integer(nk),
              K=as.integer(K),
              mu0=as.double(mu0), prior.mean.df=as.double(prior.mean.df),
              sigma0=as.double(sigma0), prior.sd.df=as.double(prior.sd.df),
              muproposal=as.double(muproposal),
              sigmaproposal=as.double(sigmaproposal),
              N=as.integer(length(pop)),
              musample=as.double(musample),
              sigmasample=as.double(sigmasample),
              samplesize=as.integer(samplesize),
              staken=integer(1),
              burnin=as.integer(burnin),
              interval=as.integer(interval),
              fVerbose=as.integer(verbose))
    Cret$sample <- cbind(Cret$musample, Cret$sigmasample)
    colnames(Cret$sample) <- c("mu","sigma")
    Cret$sample[,"mu"] <- exp(Cret$sample[,"mu"]+0.5*Cret$sample[,"sigma"]*Cret$sample[,"sigma"])
    Cret$sample[,"sigma"] <- Cret$sample[,"mu"]*sqrt(exp(Cret$sample[,"sigma"]*Cret$sample[,"sigma"])-1)
    endrun <- burnin+interval*(samplesize-1)
    attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    Cret
}
priorpln<-function(prior.mean.degree=7, prior.sd.degree=2.2,
                 prior.mean.df=4,prior.sd.df=10,
                 muproposal=0.01, 
                 sigmaproposal=0.1, 
                 samplesize=1000,burnin=0,interval=1,
                 seed=NULL,
                 verbose=TRUE){
    sigma0 <- sqrt(log(1+prior.sd.degree*prior.sd.degree/(prior.mean.degree*prior.mean.degree)))
    mu0 <- log(prior.mean.degree)-0.5*sigma0*sigma0
    musample <- rep(0,samplesize)
    sigmasample <- rep(0,samplesize)
    Cret <- .C("MHpriorpln",
              mu0=as.double(mu0), prior.mean.df=as.double(prior.mean.df),
              sigma0=as.double(sigma0), prior.sd.df=as.double(prior.sd.df),
              muproposal=as.double(muproposal),
              sigmaproposal=as.double(sigmaproposal),
              musample=as.double(musample),
              sigmasample=as.double(sigmasample),
              samplesize=as.integer(samplesize),
              staken=integer(1),
              burnin=as.integer(burnin),
              interval=as.integer(interval),
              fVerbose=as.integer(verbose))
    Cret$sample <- cbind(Cret$musample, Cret$sigmasample)
    colnames(Cret$sample) <- c("mu","sigma")
    Cret$sample[,"mu"] <- exp(Cret$sample[,"mu"]+0.5*Cret$sample[,"sigma"]*Cret$sample[,"sigma"])
    Cret$sample[,"sigma"] <- Cret$sample[,"mu"]*sqrt(exp(Cret$sample[,"sigma"]*Cret$sample[,"sigma"])-1)
    endrun <- burnin+interval*(samplesize-1)
    attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    Cret
}
posteriorsize<-function(s,
                  maxN=4*length(s),
                  K=2*max(s), nk=tabulate(s,nbin=K), n=length(s),
		  N=0.5*maxN,
                  prior.mean.degree=7, prior.sd.degree=2.2,
                  prior.mean.df=4,prior.sd.df=10,
                  muproposal=0.01, 
                  sigmaproposal=0.1, 
                  samplesize=1000,burnin=100,interval=1,burnintheta=500,
                  parallel=1, seed=NULL,
                  verbose=TRUE){
  sigma0 <- sqrt(log(1+prior.sd.degree*prior.sd.degree/(prior.mean.degree*prior.mean.degree)))
  mu0 <- log(prior.mean.degree)-0.5*sigma0*sigma0
  if(parallel==1){
      Cret <- poswest(s=s,N=N,K=K,nk=nk,n=n,maxN=maxN,
                      prior.mean.degree=prior.mean.degree,prior.mean.df=prior.mean.df,
                      prior.sd.degree=prior.sd.degree,prior.sd.df=prior.sd.df,
                      sigmaproposal=sigmaproposal, 
                      samplesize=samplesize,burnin=burnin,interval=interval,
		      burnintheta=burnintheta,
                      seed=seed)
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
    samplesize.parallel=round(samplesize/parallel)
    outlist <- clusterCall(cl, poswest,
      s=s,N=N,K=K,nk=nk,n=n,maxN=maxN,
      prior.mean.degree=prior.mean.degree,prior.mean.df=prior.mean.df,
      prior.sd.degree=prior.sd.degree,prior.sd.df=prior.sd.df,
      sigmaproposal=sigmaproposal, 
      samplesize=samplesize.parallel,burnin=burnin,interval=interval,
      burnintheta=burnintheta)
#
#   Process the results
#
    Cret <- outlist[[1]]
    for(i in (2 : length(outlist))){
     z <- outlist[[i]]
     Cret$sample <- rbind(Cret$sample,z$sample)
     Cret$samplesize <- samplesize
    }
    colnames(Cret$sample) <- c("N","mu","sigma","isolates")
    endrun <- burnin+interval*(samplesize-1)
    attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    mapfn <- function(x){
      posdensN <- density(x)
      posdensN$x[which.max(posdensN$y)]
    }
    Cret$MAP <- apply(Cret$sample,2,mapfn)
    if(verbose){
     cat("parallel samplesize=", parallel,"by", samplesize.parallel,"\n")
    }
    stopCluster(cl)
    if(getClusterOption("type")=="PVM") .PVM.exit()
  }
  Cret
}
