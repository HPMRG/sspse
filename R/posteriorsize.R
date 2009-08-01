posteriorsize<-function(s,
                  maxN=NULL,
                  K=2*max(s), n=length(s),
                  nk=tabulate(s,nbin=K),
                  mean.prior.degree=7, sd.prior.degree=3,
                  df.mean.prior=1,df.sd.prior=5,
                  muproposal=0.1, 
                  sigmaproposal=0.15, 
                  Np=0,
                  samplesize=1000,burnin=100,interval=1,burnintheta=500,
		  priorsizedistribution=c("proportion","nbinom","pln","flat"),
		  mean.prior.size=NULL, sd.prior.size=NULL,
		  mode.prior.sample.proportion=0.5,
		  median.prior.size=NULL,
		  mode.prior.size=NULL,
		  degreedistribution=c("cmp","nbinom","pln"),
                  parallel=1, seed=NULL,
                  verbose=TRUE){
#
  degreedistribution=match.arg(degreedistribution)
  posfn <- switch(degreedistribution,
                  nbinom=posnbinom,
                  pln=pospln,
		  cmp=poscmp,
		  poscmp)
  priorsizedistribution=match.arg(priorsizedistribution)
  if(priorsizedistribution=="nbinom" && missing(mean.prior.size)){
    stop("You need to specify 'mean.prior.size', and possibly 'sd.prior.size' if you use the 'nbinom' prior.") 
  }
  ### are we running the job in parallel (parallel > 1), if not just 
  #   call the degree specific function
  if(parallel==1){
      Cret <- posfn(s=s,N=N,K=K,nk=nk,n=n,maxN=maxN,
                    mean.prior.degree=mean.prior.degree,df.mean.prior=df.mean.prior,
                    sd.prior.degree=sd.prior.degree,df.sd.prior=df.sd.prior,
                    muproposal=muproposal, sigmaproposal=sigmaproposal, 
		    Np=Np,
                    samplesize=samplesize,burnin=burnin,interval=interval,
		    burnintheta=burnintheta,
		    priorsizedistribution=priorsizedistribution,
		    mean.prior.size=mean.prior.size, sd.prior.size=sd.prior.size,
		    mode.prior.sample.proportion=mode.prior.sample.proportion,
		    median.prior.size=median.prior.size,
		    mode.prior.size=mode.prior.size,
                    seed=seed)
  }
  ### since running job in parallel, start pvm (if not already running)
  else{
    cl <- beginsnow(parallel)
    ### divide the samplesize by the number of parallel runs (number of MCMC samples)
    samplesize.parallel=round(samplesize/parallel)
    ### cluster call, send following to each of the virtual machines, posnbinom function
    ### with it's arguments
    outlist <- clusterCall(cl, posfn,
      s=s,N=N,K=K,nk=nk,n=n,maxN=maxN,
      mean.prior.degree=mean.prior.degree,df.mean.prior=df.mean.prior,
      sd.prior.degree=sd.prior.degree,df.sd.prior=df.sd.prior,
      muproposal=muproposal, sigmaproposal=sigmaproposal, 
      Np=Np,
      samplesize=samplesize.parallel,burnin=burnin,interval=interval,
      burnintheta=burnintheta,
      priorsizedistribution=priorsizedistribution,
      mean.prior.size=mean.prior.size, sd.prior.size=sd.prior.size,
      mode.prior.sample.proportion=mode.prior.sample.proportion,
      median.prior.size=median.prior.size,
      mode.prior.size=mode.prior.size)
#
#   Process the results
#
    ### Snow returns a list of length parallel where each element is the return of each posfn
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
    Cret$predictive.degree<-Cret$ppos
    Cret$ppos<-NULL
    Cret$predictive.degree.count<-Cret$nk
    Cret$nk<-NULL
    #
    degnames <- NULL
    if(Np>0){degnames <- c(degnames,paste("pdeg",1:Np,sep=""))}
    colnamessample <- c("N","mu","sigma","degree1")
    if(length(degnames)>0){
     colnamessample <- c(colnamessample, degnames)
    }
    if(degreedistribution=="cmp"){
     colnamessample <- c(colnamessample, c("lambda","nu"))
    }
    colnames(Cret$sample) <- colnamessample
    
    ### Coda package which does MCMC diagnostics, requires certain attributes of MCMC sample
    endrun <- burnin+interval*(samplesize)
#!?   endrun <- burnin+interval*(samplesize-1)
    attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    
    ### define function that will compute mode of a sample
    mapfn <- function(x,lbound=min(x),ubound=max(x)){
      posdensN <- density(x, from=lbound, to=ubound)
      posdensN$x[which.max(posdensN$y)]
    }
    
    ### compute modes of posterior samples,Maximum A Posterior (MAP) values N, mu, sigma, degree1
    Cret$MAP <- apply(Cret$sample,2,mapfn)
    Cret$MAP["N"] <- mapfn(Cret$sample[,"N"],lbound=n,ubound=Cret$maxN)
    if(verbose){
     cat("parallel samplesize=", parallel,"by", samplesize.parallel,"\n")
    }
    
    ### stop cluster and PVM (in case PVM is flakey)
    endsnow(cl)
  }
  Cret$N <- c(Cret$MAP["N"], 
              mean(Cret$sample[,"N"]),
              median(Cret$sample[,"N"]),
	      quantile(Cret$sample[,"N"],c(0.025,0.975)))
  names(Cret$N) <- c("MAP","Mean AP","Median AP","P025","P975")
  #
  if(Cret$predictive.degree[length(Cret$predictive.degree)] > 0.01){
   warning("There is a non-trivial proportion of the posterior mass on very high degrees. This may indicate convergence problems in the MCMC.")
  }
  Cret$degreedistribution <- degreedistribution
  Cret$priorsizedistribution <- priorsizedistribution
  Cret$mean.prior.size <- mean.prior.size
  ### return result
  Cret
}
