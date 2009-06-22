posdis<-function(s,dis,maxN=4*length(s),
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
	      Np0i=as.integer(Np0),
	      Np1i=as.integer(Np1),
              muproposal=as.double(muproposal),
              sigmaproposal=as.double(sigmaproposal),
              N=as.integer(N),
              maxN=as.integer(maxN),
              sample=double(samplesize*dimsample),
              p0pos=double(K), p1pos=double(K),
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
    aaa <- sum(Cret$nk0+Cret$nk1)
    Cret$Nk0<-Cret$nk0/aaa
    Cret$Nk1<-Cret$nk1/aaa
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
priordis<-function(mean0.prior.degree=7,mean1.prior.degree=7,
                 sd.prior.degree=2.2,
                 df.mean.prior=1,df.sd.prior=5,
                 muproposal=0.01, 
                 sigmaproposal=0.1, 
                 samplesize=1000,burnin=0,interval=1,
                 seed=NULL,
                 verbose=TRUE){
    sigma0 <- sqrt(log(1+sd.prior.degree*sd.prior.degree/(mean0.prior.degree*mean0.prior.degree)))
    sigma1 <- sqrt(log(1+sd.prior.degree*sd.prior.degree/(mean1.prior.degree*mean1.prior.degree)))
    mu0 <- log(mean0.prior.degree)-0.5*sigma0*sigma0
    mu1 <- log(mean1.prior.degree)-0.5*sigma1*sigma1
    musample <- rep(0,2*samplesize)
    sigmasample <- rep(0,2*samplesize)
    Cret <- .C("MHpriordis",
              mu0=as.double(mu0), mu1=as.double(mu1),
	      df.mean.prior=as.double(df.mean.prior),
              sigma0=as.double(sigma0),
              sigma1=as.double(sigma1),
	      df.sd.prior=as.double(df.sd.prior),
              muproposal=as.double(muproposal),
              sigmaproposal=as.double(sigmaproposal),
              musample=as.double(musample),
              sigmasample=as.double(sigmasample),
              samplesize=as.integer(samplesize),
              staken=integer(1),
              burnin=as.integer(burnin),
              interval=as.integer(interval),
              fVerbose=as.integer(verbose))
    Cret$musample<-matrix(Cret$musample,nrow=samplesize,ncol=2,byrow=TRUE)
    Cret$sigmasample<-matrix(Cret$sigmasample,nrow=samplesize,ncol=2,byrow=TRUE)
    Cret$sample <- cbind(Cret$musample, Cret$sigmasample)
    colnames(Cret$sample) <- c("mu0","mu1","sigma0","sigma1")
    Cret$sample[,"mu0"] <- exp(Cret$sample[,"mu0"]+0.5*Cret$sample[,"sigma0"]*Cret$sample[,"sigma0"])
    Cret$sample[,"mu1"] <- exp(Cret$sample[,"mu1"]+0.5*Cret$sample[,"sigma1"]*Cret$sample[,"sigma1"])
    Cret$sample[,"sigma0"] <- Cret$sample[,"mu0"]*sqrt(exp(Cret$sample[,"sigma0"]*Cret$sample[,"sigma0"])-1)
    Cret$sample[,"sigma1"] <- Cret$sample[,"mu1"]*sqrt(exp(Cret$sample[,"sigma1"]*Cret$sample[,"sigma1"])-1)
    endrun <- burnin+interval*(samplesize-1)
    attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    Cret
}
posteriordisease<-function(s,dis,
                  maxN=4*length(s),
                  K=2*max(s), n=length(s),
		  nk0=tabulate(s[dis==0],nbin=K),
		  nk1=tabulate(s[dis==1],nbin=K),
		  N=0.5*maxN,
                  mean0.prior.degree=7,
                  mean1.prior.degree=7,
		  sd.prior.degree=2.2,
                  df.mean.prior=1,df.sd.prior=5,
                  muproposal=0.01, 
                  sigmaproposal=0.1, 
                  Np0=0, Np1=0,
                  samplesize=1000,burnin=100,interval=1,burnintheta=500,
                  parallel=1, seed=NULL,
                  verbose=TRUE){
  sigma0 <- sqrt(log(1+sd.prior.degree*sd.prior.degree/(mean0.prior.degree*mean0.prior.degree)))
  sigma1 <- sqrt(log(1+sd.prior.degree*sd.prior.degree/(mean1.prior.degree*mean1.prior.degree)))
  mu0 <- log(mean0.prior.degree)-0.5*sigma0*sigma0
  mu1 <- log(mean1.prior.degree)-0.5*sigma1*sigma1
  require(locfit, quietly=TRUE)
  mapfn <- function(x,lbound=min(x),ubound=max(x)){
    posdensN <- locfit(~ lp(x),xlim=c(lbound,ubound))
    locx <- seq(lbound,ubound,length=2000)
    locy <- predict(posdensN,newdata=locx)
    locx[which.max(locy)]
  }
  if(parallel==1){
      Cret <- posdis(s=s,dis=dis,N=N,K=K,nk0=nk0,nk1=nk1,n=n,maxN=maxN,
                      mean0.prior.degree=mean0.prior.degree,
                      mean1.prior.degree=mean1.prior.degree,
		      df.mean.prior=df.mean.prior,
                      sd.prior.degree=sd.prior.degree,df.sd.prior=df.sd.prior,
                      sigmaproposal=sigmaproposal,Np0=Np0,Np1=Np1,
                      samplesize=samplesize,burnin=burnin,interval=interval,
		      burnintheta=burnintheta,
                      seed=seed)
  }else{
    cl <- beginsnow(parallel)
    samplesize.parallel=round(samplesize/parallel)
    outlist <- clusterCall(cl, posdis,
      s=s,dis=dis,N=N,K=K,nk0=nk0,nk1=nk1,n=n,maxN=maxN,
      mean0.prior.degree=mean0.prior.degree,
      mean1.prior.degree=mean1.prior.degree,
      df.mean.prior=df.mean.prior,
      sd.prior.degree=sd.prior.degree,df.sd.prior=df.sd.prior,
      sigmaproposal=sigmaproposal,Np0=Np0,Np1=Np1,
      samplesize=samplesize.parallel,burnin=burnin,interval=interval,
      burnintheta=burnintheta)
#
#   Process the results
#
    Cret <- outlist[[1]]
    Nparallel <- length(outlist)
    Cret$samplesize <- samplesize
    for(i in (2 : Nparallel)){
     z <- outlist[[i]]
     Cret$sample <- rbind(Cret$sample,z$sample)
     Cret$Nk0<-Cret$Nk0+z$Nk0
     Cret$Nk1<-Cret$Nk1+z$Nk1
     Cret$p0pos<-Cret$p0pos+z$p0pos
     Cret$p1pos<-Cret$p1pos+z$p1pos
    }
    Cret$Nk0<-Cret$Nk0/Nparallel
    Cret$Nk1<-Cret$Nk1/Nparallel
    Cret$p0pos<-Cret$p0pos/Nparallel
    Cret$p1pos<-Cret$p1pos/Nparallel
    endrun <- burnin+interval*(samplesize)
    attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    if(verbose){
     cat("parallel samplesize=", parallel,"by", samplesize.parallel,"\n")
    }
    endsnow(cl)
  }
  degnames <- NULL
  if(Np0>0){degnames <- c(degnames,paste("p0deg",1:Np0,sep=""))}
  if(Np1>0){degnames <- c(degnames,paste("p1deg",1:Np1,sep=""))}
  colnamessample <- c("N","mu0","mu1","sigma0","sigma1","degree1",
                      "beta","disease.count")
  if(length(degnames)>0){
   colnames(Cret$sample) <- c(colnamessample,degnames,"disease")
  }else{
   colnames(Cret$sample) <- c(colnamessample,"disease")
  }
  Cret$MAP <- apply(Cret$sample,2,mapfn)
  Cret$MAP["N"] <- mapfn(Cret$sample[,"N"],lbound=n,ubound=maxN)
  Cret$N <- c(Cret$MAP["N"], 
              mean(Cret$sample[,"N"]),
              median(Cret$sample[,"N"]),
	      quantile(Cret$sample[,"N"],c(0.025,0.975)))
  Cret$disease <- c(Cret$MAP["disease"], 
              mean(Cret$sample[,"disease"]),
              median(Cret$sample[,"disease"]),
	      quantile(Cret$sample[,"disease"],c(0.025,0.975)))
  Cret
}
