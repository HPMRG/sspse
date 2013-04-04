posteriorsize<-function(s,
		  median.prior.size=NULL,
                  interval=10,
                  burnin=5000,
                  maxN=NULL,
                  K=round(quantile(s,0.80)), 
                  samplesize=1000,
		  quartiles.prior.size=NULL,
		  mean.prior.size=NULL,
		  mode.prior.size=NULL,
		  priorsizedistribution=c("proportion","flat","nbinom","pln"),
		  effective.prior.df=1,
                  sd.prior.size=NULL,
		  mode.prior.sample.proportion=0.5,
		  alpha=NULL,
		  degreedistribution=c("cmp","nbinom","pln"),
                  mean.prior.degree=NULL, sd.prior.degree=NULL,
                  df.mean.prior=1,df.sd.prior=5,
                  Np=0,
                  nk=tabulate(s,nbins=K),
                  n=length(s),
                  muproposal=0.1, 
                  sigmaproposal=0.15, 
                  burnintheta=500,
                  parallel=1, parallel.type="PVM", seed=NULL, dispersion=0,
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
  if(is.null(mean.prior.degree)){
    degs <- s
    degs[degs>K] <- K
    degs[degs==0]<-1
    isnas <- is.na(degs)
    degs <- sum(!isnas)*(degs)/sum(degs,na.rm=TRUE)
    weights <- (1/degs)
    weights[is.na(weights)] <- 0
    mean.prior.degree <- sum(s*weights)/sum(weights)
    if(is.null(sd.prior.degree)){
     sd.prior.degree <- sum(s*s*weights)/sum(weights)
     sd.prior.degree <- sqrt(sd.prior.degree - mean.prior.degree^2)
    }
  }
  cat(sprintf("The mean of the prior distribution for degree is %f.\n",mean.prior.degree))
  cat(sprintf("The s.d. of the prior distribution for degree is %f.\n",sd.prior.degree))
  if(sd.prior.degree > sqrt(4*mean.prior.degree)){
    sd.prior.degree <- min(sqrt(4*mean.prior.degree), sd.prior.degree)
    cat(sprintf("The suggested s.d. of the prior distribution for degree is too
large and has been reduced to the more reasonable %f.\n",sd.prior.degree))
  }
  ### are we running the job in parallel (parallel > 1), if not just 
  #   call the degree specific function
  if(parallel==1){
      Cret <- posfn(s=s,K=K,nk=nk,n=n,maxN=maxN,
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
		    quartiles.prior.size=quartiles.prior.size,
                    effective.prior.df=effective.prior.df,
                    alpha=alpha,
                    seed=seed,
                    dispersion=dispersion)
  }
  ### since running job in parallel, start pvm (if not already running)
  else{
    cl <- beginparallel(parallel,type=parallel.type)
    ### divide the samplesize by the number of parallel runs (number of MCMC samples)
    samplesize.parallel=round(samplesize/parallel)
    ### cluster call, send following to each of the virtual machines, posnbinom function
    ### with it's arguments
    outlist <- parallel::clusterCall(cl, posfn,
      s=s,K=K,nk=nk,n=n,maxN=maxN,
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
      mode.prior.size=mode.prior.size,
      quartiles.prior.size=quartiles.prior.size,
      effective.prior.df=effective.prior.df,
      alpha=alpha,
      dispersion=dispersion)
#
#   Process the results
#
    ### Snow returns a list of length parallel where each element is the return of each posfn
    ### Following loops combines the separate MCMC samples into 1 using rbind, creating a matrix
    Cret <- outlist[[1]]
    Nparallel <- length(outlist)
    Cret$samplesize <- samplesize
    for(i in (2 : Nparallel)){
     z <- outlist[[i]]
     Cret$sample <- rbind(Cret$sample,z$sample)
     Cret$predictive.degree.count<-Cret$predictive.degree.count+z$predictive.degree.count
     Cret$predictive.degree<-Cret$predictive.degree+z$predictive.degree
    }
    Cret$predictive.degree.count<-Cret$predictive.degree.count/Nparallel
    Cret$predictive.degree<-Cret$predictive.degree/Nparallel
    #
    degnames <- NULL
    if(Np>0){degnames <- c(degnames,paste("pdeg",1:Np,sep=""))}
    colnamessample <- c("N","mu","sigma","degree1","totalsize")
    if(length(degnames)>0){
     colnamessample <- c(colnamessample, degnames)
    }
    if(degreedistribution=="cmp"){
     colnamessample <- c(colnamessample, c("lambda","nu"))
    }
    colnames(Cret$sample) <- colnamessample
    
    ### Coda package which does MCMC diagnostics, requires certain attributes of MCMC sample
    endrun <- burnin+interval*(samplesize)
    attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    
#   ### Remove the padding from the last draws from the populations of degrees
#   Nlastpos <- Cret$sample[nrow(Cret$sample),"N"]
#   Cret$pop<-Cret$pop[1:Nlastpop]

    ### compute modes of posterior samples,Maximum A Posterior (MAP) values N, mu, sigma, degree1
    Cret$MAP <- apply(Cret$sample,2,mode.density)
    Cret$MAP["N"] <- mode.density(Cret$sample[,"N"],lbound=n,ubound=Cret$maxN)
    if(verbose){
     cat("parallel samplesize=", parallel,"by", samplesize.parallel,"\n")
    }
    
    ### stop cluster and PVM (in case PVM is flakey)
    endparallel(cl,type=parallel.type)
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
# Cret$mean.prior.size <- mean.prior.size
  ### return result
  Cret
}
