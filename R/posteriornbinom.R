posteriornbinom<-function(s,
                  maxN=4*length(s),
                  K=2*max(s), nk=tabulate(s,nbin=K), n=length(s),
		  N=0.5*maxN,
                  mean.prior.degree=7, sd.prior.degree=3,
                  df.mean.prior=1,df.sd.prior=5,
                  muproposal=0.1, 
                  sigmaproposal=0.15, 
                  Np=0,
                  samplesize=1000,burnin=100,interval=1,burnintheta=500,
		  mean.prior.size=N, sd.prior.size=ceiling(sqrt(mean.prior.size)*3),
                  parallel=1, seed=NULL,
                  verbose=TRUE){
  
  ### are we running the job in parallel (parallel > 1), if not just call posnbinom
  if(parallel==1){
      Cret <- posnbinom(s=s,N=N,K=K,nk=nk,n=n,maxN=maxN,
                      mean.prior.degree=mean.prior.degree,df.mean.prior=df.mean.prior,
                      sd.prior.degree=sd.prior.degree,df.sd.prior=df.sd.prior,
                      sigmaproposal=sigmaproposal, Np=Np,
                      samplesize=samplesize,burnin=burnin,interval=interval,
		      burnintheta=burnintheta,
		      mean.prior.size=mean.prior.size, sd.prior.size=sd.prior.size,
                      seed=seed)
  }
  ### since running job in parallel, start pvm (if not already running)
  else{
    cl <- beginsnow(parallel)
    ### divide the samplesize by the number of parallel runs (number of MCMC samples)
    samplesize.parallel=round(samplesize/parallel)
    ### cluster call, send following to each of the virtual machines, posnbinom function
    ### with it's arguments
    outlist <- clusterCall(cl, posnbinom,
      s=s,N=N,K=K,nk=nk,n=n,maxN=maxN,
      mean.prior.degree=mean.prior.degree,df.mean.prior=df.mean.prior,
      sd.prior.degree=sd.prior.degree,df.sd.prior=df.sd.prior,
      sigmaproposal=sigmaproposal, Np=Np,
      samplesize=samplesize.parallel,burnin=burnin,interval=interval,
      burnintheta=burnintheta,
      mean.prior.size=mean.prior.size, sd.prior.size=sd.prior.size)
#
#   Process the results
#
    ### Snow returns a list of length parallel where each element is the return of each posnbinom
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
    mapfn <- function(x,lbound=min(x),ubound=max(x)){
      posdensN <- density(x, from=lbound, to=ubound)
      posdensN$x[which.max(posdensN$y)]
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
posnbinom<-function(s,maxN=4*length(s),
                  K=2*max(s), nk=tabulate(s,nbin=K), n=length(s),
		  N=maxN/2,
                  mean.prior.degree=7, sd.prior.degree=3,
                  df.mean.prior=1, df.sd.prior=5,
                  muproposal=0.25, 
                  sigmaproposal=0.15, 
                  Np=0,
                  samplesize=10,burnin=0,interval=1,burnintheta=500,
		  mean.prior.size=N, sd.prior.size=ceiling(sqrt(mean.prior.size)*3),
                  seed=NULL,
                  verbose=TRUE){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    sigma0 <- sd.prior.degree
    mu0 <- mean.prior.degree
    if(!is.null(seed))  set.seed(as.integer(seed))
    dimsample <- 4+Np
    lpriorm <- dnbinommu(x=n+(1:maxN)-1,
                         mu=mean.prior.size, sd=sd.prior.size,
			 log=TRUE)
 #print(cbind(n+(1:maxN)-1,lpriorm))
    Cret <- .C("gnbinom",
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
              lpriorm=as.double(lpriorm),
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
    Cret$Nk<-Cret$nk/sum(Cret$nk)
    endrun <- burnin+interval*(samplesize-1)
    attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    require(locfit, quietly=TRUE)
    mapfn <- function(x,lbound=min(x),ubound=max(x)){
      posdensN <- locfit(~ lp(x),xlim=c(lbound,ubound))
      locx <- seq(lbound,ubound,length=2000)
      locy <- predict(posdensN,newdata=locx)
      locx[which.max(locy)]
    }
    mapfn <- function(x,lbound=min(x),ubound=max(x)){
      posdensN <- density(x, from=lbound, to=ubound)
      posdensN$x[which.max(posdensN$y)]
    }
    Cret$MAP <- apply(Cret$sample,2,mapfn)
    Cret$MAP["N"] <- mapfn(Cret$sample[,"N"],lbound=n,ubound=maxN)
    Cret
}

dnbinommu <- function(x, mu, sd, log=FALSE){
    if (missing(mu)) {
      stop("The mean of the Negative Binomial must be specified.")
    }
    if (missing(sd)) {
      stop("The standard deviation of the Negative Binomial must be specified.")
    }
   out <- .C("dnbmu",
            x=as.integer(x),
            mu=as.double(mu),
            sig=as.double(sd),
            n=as.integer(length(x)),
            give_log=as.integer(log),
            val=double(length(x)),
            PACKAGE="size")
   return(out$val)
}
