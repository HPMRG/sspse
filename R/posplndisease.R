posplndisease<-function(s,dis,maxN=4*length(s),
                  K=2*max(s), 
		  nk0=tabulate(s[dis==0],nbin=K),
		  nk1=tabulate(s[dis==1],nbin=K),
		  n=length(s),
		  N=maxN/2,
                  mean0.prior.degree=7, 
                  mean1.prior.degree=7, 
		  sd.prior.degree=3,
                  df.mean.prior=1, df.sd.prior=5,
                  muproposal=0.1, 
                  sigmaproposal=0.15, 
                  Np0=0, Np1=0,
                  samplesize=10,burnin=0,interval=1,burnintheta=500,
		  mean.prior.size=N, sd.prior.size=N,
                  seed=NULL,
                  verbose=TRUE){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    if(!is.null(seed))  set.seed(as.integer(seed))
    #
    # Transform observed mean parametrization to log-normal
    # parametrization
    #
    sigma0 <- sqrt(log(1+sd.prior.degree*sd.prior.degree/(mean0.prior.degree*mean0.prior.degree)))
    sigma1 <- sqrt(log(1+sd.prior.degree*sd.prior.degree/(mean1.prior.degree*mean1.prior.degree)))
    mu0 <- log(mean0.prior.degree)-0.5*sigma0*sigma0
    mu1 <- log(mean1.prior.degree)-0.5*sigma1*sigma1
    dimsample <- 8+Np0+Np1
    lpriorm <- dnbinommu(x=n+(1:maxN)-1,
                         mu=mean.prior.size, sd=sd.prior.size,
			 log=TRUE)
    Cret <- .C("gplndisease",
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
              lpriorm=as.double(lpriorm),
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
    #
    # Transform observed mean parametrization to log-normal
    # parametrization
    #
    # Expectation and s.d. of normal from log-normal
    #
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
    mapfn <- function(x,lbound=min(x),ubound=max(x)){
      posdensN <- density(x, from=lbound, to=ubound)
      posdensN$x[which.max(posdensN$y)]
    }
    Cret$MAP <- apply(Cret$sample,2,mapfn)
    Cret$MAP["N"] <- mapfn(Cret$sample[,"N"],lbound=n,ubound=maxN)
    Cret
}
