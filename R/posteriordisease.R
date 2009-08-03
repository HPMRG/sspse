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
		  degreedistribution=c("cmp","nbinom","pln"),
                  maxN=NULL,
                  K=2*max(s), n=length(s),
		  nk0=tabulate(s[dis==0],nbin=K),
		  nk1=tabulate(s[dis==1],nbin=K),
                  muproposal=0.1, 
                  sigmaproposal=0.15, 
                  parallel=1, seed=NULL,
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
                      seed=seed)
  }else{
    cl <- beginsnow(parallel)
    samplesize.parallel=round(samplesize/parallel)
    outlist <- clusterCall(cl, posfn,
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
      mode.prior.size=mode.prior.size)
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
posdis.old<-function(s,dis,maxN=4*length(s),
                  K=2*max(s), 
		  nk0=tabulate(s[dis==0],nbin=K),
		  nk1=tabulate(s[dis==1],nbin=K),
		  n=length(s),
		  N=maxN/2,
                  mean0.prior.degree=7, 
                  mean1.prior.degree=7, 
		  sd.prior.degree=2.2,
                  df.mean.prior=1, df.sd.prior=5,
                  muproposal=0.01, 
                  sigmaproposal=0.1, 
                  Np0=0, Np1=0,
                  samplesize=10,burnin=0,interval=1,burnintheta=500,
                  seed=NULL,
                  verbose=TRUE){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    sigma0 <- sqrt(log(1+sd.prior.degree*sd.prior.degree/(mean0.prior.degree*mean0.prior.degree)))
    sigma1 <- sqrt(log(1+sd.prior.degree*sd.prior.degree/(mean1.prior.degree*mean1.prior.degree)))
    mu0 <- log(mean0.prior.degree)-0.5*sigma0*sigma0
    mu1 <- log(mean1.prior.degree)-0.5*sigma1*sigma1
    if(!is.null(seed))  set.seed(as.integer(seed))
    dimsample <- 8+Np0+Np1
    Cret <- .C("gsppsdis",
              pop=as.integer(c(s,rep(0,maxN-n))),
              dis=as.integer(c(dis,rep(0,maxN-n))),
              nk0=as.integer(nk0),
              nk1=as.integer(nk1),
              K=as.integer(K),
              n=as.integer(n),
              samplesize=as.integer(samplesize),
              burnin=as.integer(burnin),
              interval=as.integer(interval),
              mu0=as.double(mu0), mu1=as.double(mu1),
	      df.mean.prior=as.double(df.mean.prior),
              sigma0=as.double(sigma0),
              sigma1=as.double(sigma1),
	      df.sd.prior=as.double(df.sd.prior),
	      Np0=as.integer(Np0),
	      Np1=as.integer(Np1),
              muproposal=as.double(muproposal),
              sigmaproposal=as.double(sigmaproposal),
              N=as.integer(N),
              maxN=as.integer(maxN),
              sample=double(samplesize*dimsample),
              p0pos=double(K), p1pos=double(K),
              ppos=double(K),
              burnintheta=as.integer(burnintheta),
              fVerbose=as.integer(verbose))
    Cret$sample<-matrix(Cret$sample,nrow=samplesize,ncol=dimsample,byrow=TRUE)
    degnames <- NULL
    if(Np0>0){degnames <- c(degnames,paste("p0deg",1:Np0,sep=""))}
    if(Np1>0){degnames <- c(degnames,paste("p1deg",1:Np1,sep=""))}
    colnamessample <- c("N","mu0","mu1","sigma0","sigma1","degree1",
                        "beta","disease.count")
    if(length(degnames)>0){
     colnames(Cret$sample) <- c(colnamessample, degnames)
    }else{
     colnames(Cret$sample) <- colnamessample
    }
    Cret$sample<-cbind(Cret$sample,Cret$sample[,"disease.count"]/Cret$sample[,"N"])
    if(length(degnames)>0){
     colnames(Cret$sample) <- c(colnamessample,degnames,"disease")
    }else{
     colnames(Cret$sample) <- c(colnamessample,"disease")
    }
    # Expectation and s.d. of log-normal
    Cret$sample[,"mu0"] <- exp(Cret$sample[,"mu0"]+0.5*Cret$sample[,"sigma0"]*Cret$sample[,"sigma0"])
    Cret$sample[,"mu1"] <- exp(Cret$sample[,"mu1"]+0.5*Cret$sample[,"sigma1"]*Cret$sample[,"sigma1"])
    Cret$sample[,"sigma0"] <- Cret$sample[,"mu0"]*sqrt(exp(Cret$sample[,"sigma0"]*Cret$sample[,"sigma0"])-1)
    Cret$sample[,"sigma1"] <- Cret$sample[,"mu1"]*sqrt(exp(Cret$sample[,"sigma1"]*Cret$sample[,"sigma1"])-1)
    # Expectation and s.d. of Poisson-log-normal
    Cret$sample[,"sigma0"] <- sqrt(Cret$sample[,"mu0"]+Cret$sample[,"sigma0"]*Cret$sample[,"sigma0"])
    Cret$sample[,"sigma1"] <- sqrt(Cret$sample[,"mu1"]+Cret$sample[,"sigma1"]*Cret$sample[,"sigma1"])
#   aaa <- sum(Cret$nk0+Cret$nk1)
#   Cret$Nk0<-Cret$nk0/aaa
#   Cret$Nk1<-Cret$nk1/aaa
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
    Cret$MAP <- apply(Cret$sample,2,mapfn)
    Cret$MAP["N"] <- mapfn(Cret$sample[,"N"],lbound=n,ubound=maxN)
    Cret
}
