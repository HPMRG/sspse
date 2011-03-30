posteriordisease<-function(s,dis,
                  mean0.prior.degree=7,
                  mean1.prior.degree=7,
		  sd.prior.degree=3,
                  df.mean.prior=1,df.sd.prior=5,
                  Np0=0, Np1=0,
                  samplesize=1000,burnin=100,interval=1,burnintheta=500,
		  priorsizedistribution=c("proportion","nbinom","pln","flat"),
		  mean.prior.size=NULL, sd.prior.size=NULL,
		  mode.prior.sample.proportion=0.5,
		  median.prior.size=NULL,
		  mode.prior.size=NULL,
		  effective.prior.df=1,
		  degreedistribution=c("cmp","nbinom","pln"),
                  maxN=NULL,
                  K=round(quantile(s,0.80)), n=length(s),
		  nk0=tabulate(s[dis==0],nbin=K),
		  nk1=tabulate(s[dis==1],nbin=K),
                  muproposal=0.1, 
                  sigmaproposal=0.15, 
                  parallel=1, seed=NULL, dispersion=0,
                  verbose=TRUE){
  degreedistribution=match.arg(degreedistribution)
  posfn <- switch(degreedistribution,
                  nbinom=posnbinomdisease,
                  pln=posplndisease,
		  cmp=poscmpdisease,
		  poscmpdisease)
  priorsizedistribution=match.arg(priorsizedistribution)
  if(priorsizedistribution=="nbinom" && missing(mean.prior.size)){
    stop("You need to specify 'mean.prior.size', and possibly 'sd.prior.size' if you use the 'nbinom' prior.") 
  }
  cat("K = ",K,"\n")
  if(parallel==1){
      Cret <- posfn(s=s,dis=dis,K=K,nk0=nk0,nk1=nk1,n=n,maxN=maxN,
                      mean0.prior.degree=mean0.prior.degree,
                      mean1.prior.degree=mean1.prior.degree,
		      df.mean.prior=df.mean.prior,
                      sd.prior.degree=sd.prior.degree,df.sd.prior=df.sd.prior,
		      muproposal=muproposal, sigmaproposal=sigmaproposal,
                      Np0=Np0,Np1=Np1,
                      samplesize=samplesize,burnin=burnin,interval=interval,
		      burnintheta=burnintheta,
		      priorsizedistribution=priorsizedistribution,
		      mean.prior.size=mean.prior.size, sd.prior.size=sd.prior.size,
		      mode.prior.sample.proportion=mode.prior.sample.proportion,
		      median.prior.size=median.prior.size,
		      mode.prior.size=mode.prior.size,
                      effective.prior.df=effective.prior.df,
                      seed=seed,
		      dispersion=dispersion)
  }else{
    cl <- beginsnow(parallel)
    samplesize.parallel=round(samplesize/parallel)
    outlist <- snow::clusterCall(cl, posfn,
      s=s,dis=dis,K=K,nk0=nk0,nk1=nk1,n=n,maxN=maxN,
      mean0.prior.degree=mean0.prior.degree,
      mean1.prior.degree=mean1.prior.degree,
      df.mean.prior=df.mean.prior,
      sd.prior.degree=sd.prior.degree,df.sd.prior=df.sd.prior,
      muproposal=muproposal, sigmaproposal=sigmaproposal,
      Np0=Np0,Np1=Np1,
      samplesize=samplesize.parallel,burnin=burnin,interval=interval,
      burnintheta=burnintheta,
      priorsizedistribution=priorsizedistribution,
      mean.prior.size=mean.prior.size, sd.prior.size=sd.prior.size,
      mode.prior.sample.proportion=mode.prior.sample.proportion,
      median.prior.size=median.prior.size,
      mode.prior.size=mode.prior.size,
      effective.prior.df=effective.prior.df,
      dispersion=dispersion)
#
#   Process the results
#
    Cret <- outlist[[1]]
    Nparallel <- length(outlist)
    Cret$samplesize <- samplesize
    for(i in (2 : Nparallel)){
     z <- outlist[[i]]
     Cret$sample <- rbind(Cret$sample,z$sample)
     Cret$predictive.degree.count0<-Cret$predictive.degree.count0+z$predictive.degree.count0
     Cret$predictive.degree.count1<-Cret$predictive.degree.count1+z$predictive.degree.count1
     Cret$predictive.degree0<-Cret$predictive.degree0+z$predictive.degree0
     Cret$predictive.degree1<-Cret$predictive.degree1+z$predictive.degree1
     Cret$predictive.degree<-Cret$predictive.degree+z$predictive.degree
    }
    Cret$predictive.degree.count0<-Cret$predictive.degree.count0/Nparallel
    Cret$predictive.degree.count1<-Cret$predictive.degree.count1/Nparallel
    Cret$predictive.degree0<-Cret$predictive.degree0/Nparallel
    Cret$predictive.degree1<-Cret$predictive.degree1/Nparallel
    Cret$predictive.degree<-Cret$predictive.degree/Nparallel
    #
    degnames <- NULL
    if(Np0>0){degnames <- c(degnames,paste("p0deg",1:Np0,sep=""))}
    if(Np1>0){degnames <- c(degnames,paste("p1deg",1:Np1,sep=""))}
    colnamessample <- c("N","mu0","mu1","sigma0","sigma1","degree1",
                        "beta","disease.count")
    if(length(degnames)>0){
     colnamessample <- c(colnamessample,degnames,"disease")
    }else{
     colnamessample <- c(colnamessample,"disease")
    }
    if(degreedistribution=="cmp"){
     colnamessample <- c(colnamessample, c("lambda0","nu0","lambda1","nu1"))
    }
    colnames(Cret$sample) <- colnamessample
    endrun <- burnin+interval*(samplesize)
    attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    ### define function that will compute mode of a sample
    mapfn <- function(x,lbound=min(x),ubound=max(x)){
      posdensN <- density(x, from=lbound, to=ubound)
      posdensN$x[which.max(posdensN$y)]
    }
    Cret$MAP <- apply(Cret$sample,2,mapfn)
    Cret$MAP["N"] <- mapfn(Cret$sample[,"N"],lbound=n,ubound=Cret$maxN)
    if(verbose){
     cat("parallel samplesize=", parallel,"by", samplesize.parallel,"\n")
    }
    endsnow(cl)
  }
  Cret$N <- c(Cret$MAP["N"], 
              mean(Cret$sample[,"N"]),
              median(Cret$sample[,"N"]),
	      quantile(Cret$sample[,"N"],c(0.025,0.975)))
  names(Cret$N) <- c("MAP","Mean AP","Median AP","P025","P975")
  Cret$disease <- c(Cret$MAP["disease"], 
              mean(Cret$sample[,"disease"]),
              median(Cret$sample[,"disease"]),
	      quantile(Cret$sample[,"disease"],c(0.025,0.975)))
  names(Cret$disease) <- c("MAP","Mean AP","Median AP","P025","P975")
  Cret$disease.count <- c(Cret$MAP["disease.count"], 
              mean(Cret$sample[,"disease.count"]),
              median(Cret$sample[,"disease.count"]),
	      quantile(Cret$sample[,"disease.count"],c(0.025,0.975)))
  names(Cret$disease.count) <- c("MAP","Mean AP","Median AP","P025","P975")
  #
  if(Cret$predictive.degree0[length(Cret$predictive.degree0)] > 0.01){
   warning("There is a non-trivial proportion of the posterior mass on very high degrees for non-diseased nodes. This may indicate convergence problems in the MCMC.")
  }
  if(Cret$predictive.degree1[length(Cret$predictive.degree1)] > 0.01){
   warning("There is a non-trivial proportion of the posterior mass on very high degrees for diseased nodes. This may indicate convergence problems in the MCMC.")
  }
  if(Cret$predictive.degree[length(Cret$predictive.degree)] > 0.01){
   warning("There is a non-trivial proportion of the posterior mass on very high degree nodes. This may indicate convergence problems in the MCMC.")
  }
  Cret$degreedistribution <- degreedistribution
  Cret$priorsizedistribution <- priorsizedistribution
  Cret$mean.prior.size <- mean.prior.size
  ### return result
  Cret
}
