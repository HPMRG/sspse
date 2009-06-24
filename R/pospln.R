posteriorsize<-function(s,
                  maxN=4*length(s),
                  K=2*max(s), nk=tabulate(s,nbin=K), n=length(s),
		  N=0.5*maxN,
                  mean.prior.degree=7, sd.prior.degree=3,
                  df.mean.prior=1,df.sd.prior=5,
                  muproposal=0.25, 
                  sigmaproposal=0.15, 
                  Np=0,
                  samplesize=1000,burnin=100,interval=1,burnintheta=500,
                  parallel=1, seed=NULL,
                  verbose=TRUE){
  ### takes mean and standard deviation of the prior lognormal distribution and computes the corresponding
  ### mean and standard deviation of the underlying normal
  # sigma0 <- sqrt(log(1+(sd.prior.degree*sd.prior.degree-mean.prior.degree)/(mean.prior.degree*mean.prior.degree)))
  # mu0 <- log(mean.prior.degree)-0.5*sigma0*sigma0
  
  ### are we running the job in parallel (parallel > 1), if not just call poswest
  if(parallel==1){
      Cret <- poswest(s=s,N=N,K=K,nk=nk,n=n,maxN=maxN,
                      mean.prior.degree=mean.prior.degree,df.mean.prior=df.mean.prior,
                      sd.prior.degree=sd.prior.degree,df.sd.prior=df.sd.prior,
                      sigmaproposal=sigmaproposal, Np=Np,
                      samplesize=samplesize,burnin=burnin,interval=interval,
		      burnintheta=burnintheta,
                      seed=seed)
  }
  ### since running job in parallel, start pvm (if not already running)
  else{
    cl <- beginsnow(parallel)
    ### divide the samplesize by the number of parallel runs (number of MCMC samples)
    samplesize.parallel=round(samplesize/parallel)
    ### cluster call, send following to each of the virtual machines, poswest function
    ### with it's arguments
    outlist <- clusterCall(cl, poswest,
      s=s,N=N,K=K,nk=nk,n=n,maxN=maxN,
      mean.prior.degree=mean.prior.degree,df.mean.prior=df.mean.prior,
      sd.prior.degree=sd.prior.degree,df.sd.prior=df.sd.prior,
      sigmaproposal=sigmaproposal, Np=Np,
      samplesize=samplesize.parallel,burnin=burnin,interval=interval,
      burnintheta=burnintheta)
#
#   Process the results
#
    ### Snow returns a list of length parallel where each element is the return of each poswest
    ### Following loops combines the separate MCMC samples into 1 using rbind, creating a matrix
    Cret <- outlist[[1]]
    Cret$samplesize <- samplesize
    Nparallel <- length(outlist)
    for(i in (2 : Nparallel)){
     z <- outlist[[i]]
     Cret$sample <- rbind(Cret$sample,z$sample)
     Cret$nk<-Cret$nk+z$nk
     Cret$ppos<-Cret$ppos+z$ppos
    }
    Cret$nk<-Cret$nk/Nparallel
    Cret$ppos<-Cret$ppos/Nparallel
    degnames <- NULL
    if(Np>0){degnames <- c(degnames,paste("pdeg",1:Np,sep=""))}
    colnamessample <- c("N","mu","sigma","degree1")
    if(length(degnames)>0){
     colnames(Cret$sample) <- c(colnamessample, degnames)
    }else{
     colnames(Cret$sample) <- colnamessample
    }
    
    ### Coda package which does MCMC diagnostics, requires certain attributes of MCMC sample
    endrun <- burnin+interval*(samplesize-1)
    attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    
    ### define function that will compute mode of a sample
    require(locfit, quietly=TRUE)
    mapfn <- function(x,lbound=min(x),ubound=max(x)){
      posdensN <- locfit(~ lp(x,maxk=500),xlim=c(lbound,ubound))
      locx <- seq(lbound,ubound,length=2000)
      locy <- predict(posdensN,newdata=locx)
      locx[which.max(locy)]
    }
    
    ### compute modes of posterior samples,Maximum A Posterior (MAP) values N, mu, sigma, degree1
    Cret$MAP <- apply(Cret$sample,2,mapfn)
    Cret$MAP["N"] <- mapfn(Cret$sample[,"N"],lbound=n,ubound=maxN)
    if(verbose){
     cat("parallel samplesize=", parallel,"by", samplesize.parallel,"\n")
    }
    
    ### stop cluster and PVM (in case PVM is flakey)
    endsnow(cl)
  }
  if(Cret$ppos[length(Cret$ppos)] > 0.01){
   warning("There is a non-trivial proportion of the posterior mass on very high degrees. This may indicate convergence problems in the MCMC.")
  }
  ### return result
  Cret
}
poswest<-function(s,maxN=4*length(s),
                  K=2*max(s), nk=tabulate(s,nbin=K), n=length(s),
		  N=maxN/2,
                  mean.prior.degree=7, sd.prior.degree=3,
                  df.mean.prior=1, df.sd.prior=5,
                  muproposal=0.25, 
                  sigmaproposal=0.15, 
                  Np=0,
                  samplesize=10,burnin=0,interval=1,burnintheta=500,
                  seed=NULL,
                  verbose=TRUE){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    sigma0 <- sqrt(log(1+(sd.prior.degree*sd.prior.degree-mean.prior.degree)/(mean.prior.degree*mean.prior.degree)))
    mu0 <- log(mean.prior.degree)-0.5*sigma0*sigma0
    if(!is.null(seed))  set.seed(as.integer(seed))
    dimsample <- 4+Np
    Cret <- .C("gspps",
              pop=as.integer(c(s,rep(0,maxN-n))),
              nk=as.integer(nk),
              K=as.integer(K),
              n=as.integer(n),
              samplesize=as.integer(samplesize),
              burnin=as.integer(burnin),
              interval=as.integer(interval),
              mu0=as.double(mu0), df.mean.prior=as.double(df.mean.prior),
              sigma0=as.double(sigma0), df.sd.prior=as.double(df.sd.prior),
              Npi=as.integer(Np),
              muproposal=as.double(muproposal),
              sigmaproposal=as.double(sigmaproposal),
              N=as.integer(N),
              maxN=as.integer(maxN),
              sample=double(samplesize*dimsample),
              ppos=double(K),
              burnintheta=as.integer(burnintheta),
              fVerbose=as.integer(verbose))
    Cret$sample<-matrix(Cret$sample,nrow=samplesize,ncol=dimsample,byrow=TRUE)
    degnames <- NULL
    if(Np>0){degnames <- c(degnames,paste("pdeg",1:Np,sep=""))}
    colnamessample <- c("N","mu","sigma","degree1")
    if(length(degnames)>0){
     colnames(Cret$sample) <- c(colnamessample, degnames)
    }else{
     colnames(Cret$sample) <- colnamessample
    }
    # Expectation and s.d. of normal from log-normal
    Cret$sample[,"mu"] <- exp(Cret$sample[,"mu"]+0.5*Cret$sample[,"sigma"]*Cret$sample[,"sigma"])
    Cret$sample[,"sigma"] <- Cret$sample[,"mu"]*sqrt(exp(Cret$sample[,"sigma"]*Cret$sample[,"sigma"])-1)
    # Expectation and s.d. of Poisson-log-normal
    Cret$sample[,"sigma"] <- sqrt(Cret$sample[,"mu"]+Cret$sample[,"sigma"]*Cret$sample[,"sigma"])
    Cret$Nk<-Cret$nk/sum(Cret$nk)
    endrun <- burnin+interval*(samplesize-1)
    attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    mapfn <- function(x){
      posdensN <- density(x)
      posdensN$x[which.max(posdensN$y)]
    }
    require(locfit, quietly=TRUE)
    mapfn <- function(x,lbound=min(x),ubound=max(x)){
      posdensN <- locfit(~ lp(x),xlim=c(lbound,ubound))
      locx <- seq(lbound,ubound,length=2000)
      locy <- predict(posdensN,newdata=locx)
      locx[which.max(locy)]
    }
    Cret$MAP <- apply(Cret$sample,2,mapfn)
    Cret$MAP["N"] <- mapfn(Cret$sample[,"N"],lbound=n,ubound=maxN)
#   Cret$MAP <- apply(Cret$sample,2,mapfn,lbound=n,ubound=maxN)
    Cret
}
poswestN<-function(s,N,
                 K=max(s), nk=tabulate(s,nbin=K), n=length(s),
                 mean.prior.degree=7, sd.prior.degree=3,
                 df.mean.prior=1, df.sd.prior=5,
                 muproposal=0.25, 
                 sigmaproposal=0.15, 
                 samplesize=10,burnin=0,interval=1,burnintheta=500,
                 seed=NULL,
                 verbose=TRUE){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    sigma0 <- sqrt(log(1+(sd.prior.degree*sd.prior.degree-mean.prior.degree)/(mean.prior.degree*mean.prior.degree)))
    mu0 <- log(mean.prior.degree)-0.5*sigma0*sigma0
    if(!is.null(seed))  set.seed(as.integer(seed))
    Cret <- .C("gsppsN",
              pop=as.integer(c(s,rep(0,N-n))),
              nk=as.integer(nk),
              K=as.integer(K),
              n=as.integer(n),
              samplesize=as.integer(samplesize),
              burnin=as.integer(burnin),
              interval=as.integer(interval),
              mu0=as.double(mu0), df.mean.prior=as.double(df.mean.prior),
              sigma0=as.double(sigma0), df.sd.prior=as.double(df.sd.prior),
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
                 mean.prior.degree=7, sd.prior.degree=3,
                 df.mean.prior=1, df.sd.prior=5,
                 muproposal=0.25, 
                 sigmaproposal=0.15, 
                 samplesize=1000,burnin=0,interval=1,
                 seed=NULL,
                 verbose=TRUE){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    sigma0 <- sqrt(log(1+(sd.prior.degree*sd.prior.degree-mean.prior.degree)/(mean.prior.degree*mean.prior.degree)))
    mu0 <- log(mean.prior.degree)-0.5*sigma0*sigma0
    if(!is.null(seed))  set.seed(as.integer(seed))
    Nk <- tabulate(pop,nbin=K)
    musample <- rep(0,samplesize)
    sigmasample <- rep(0,samplesize)
    Cret <- .C("MHpln",
              Nk=as.integer(Nk),
              K=as.integer(K),
              mu0=as.double(mu0), df.mean.prior=as.double(df.mean.prior),
              sigma0=as.double(sigma0), df.sd.prior=as.double(df.sd.prior),
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
    # Expectation and s.d. of normal from log-normal
    Cret$sample[,"mu"] <- exp(Cret$sample[,"mu"]+0.5*Cret$sample[,"sigma"]*Cret$sample[,"sigma"])
    Cret$sample[,"sigma"] <- Cret$sample[,"mu"]*sqrt(exp(Cret$sample[,"sigma"]*Cret$sample[,"sigma"])-1)
    # Expectation and s.d. of Poisson-log-normal
    Cret$sample[,"sigma"] <- sqrt(Cret$sample[,"mu"]+Cret$sample[,"sigma"]*Cret$sample[,"sigma"])
    endrun <- burnin+interval*(samplesize-1)
    attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    Cret
}
priorpln<-function(mean.prior.degree=7, sd.prior.degree=3,
                 df.mean.prior=1,df.sd.prior=5,
                 muproposal=0.25, 
                 sigmaproposal=0.15, 
                 samplesize=1000,burnin=0,interval=1,
                 seed=NULL,
                 verbose=TRUE){
    sigma0 <- sqrt(log(1+(sd.prior.degree*sd.prior.degree-mean.prior.degree)/(mean.prior.degree*mean.prior.degree)))
    mu0 <- log(mean.prior.degree)-0.5*sigma0*sigma0
    musample <- rep(0,samplesize)
    sigmasample <- rep(0,samplesize)
    Cret <- .C("MHpriorpln",
              mu0=as.double(mu0), df.mean.prior=as.double(df.mean.prior),
              sigma0=as.double(sigma0), df.sd.prior=as.double(df.sd.prior),
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
    # Expectation and s.d. of log-normal
    Cret$sample[,"mu"] <- exp(Cret$sample[,"mu"]+0.5*Cret$sample[,"sigma"]*Cret$sample[,"sigma"])
    Cret$sample[,"sigma"] <- Cret$sample[,"mu"]*sqrt(exp(Cret$sample[,"sigma"]*Cret$sample[,"sigma"])-1)
    # Expectation and s.d. of Poisson-log-normal
    Cret$sample[,"sigma"] <- sqrt(Cret$sample[,"mu"]+Cret$sample[,"sigma"]*Cret$sample[,"sigma"])
    endrun <- burnin+interval*(samplesize-1)
    attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    Cret
}
