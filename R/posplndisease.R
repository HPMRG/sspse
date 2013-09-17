posplndisease<-function(s,dis,
                  maxN=NULL,
                  K=2*max(s), 
		  nk0=tabulate(s[dis==0],nbins=K),
		  nk1=tabulate(s[dis==1],nbins=K),
		  n=length(s),
                  mean0.prior.degree=7, 
                  mean1.prior.degree=7, 
		  sd.prior.degree=3,
                  df.mean.prior=1, df.sd.prior=5,
                  muproposal=0.1, 
                  sigmaproposal=0.15, 
                  Np0=0, Np1=0,
                  samplesize=10,burnin=0,interval=1,burnintheta=500,
		  priorsizedistribution=c("beta","nbinom","pln","flat"),
		  mean.prior.size=NULL, sd.prior.size=NULL,
		  mode.prior.sample.proportion=0.5,
		  median.prior.sample.proportion=NULL,
		  median.prior.size=NULL,
		  mode.prior.size=NULL,
		  quartiles.prior.size=NULL,
		  effective.prior.df=1,
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
    #
    priorsizedistribution=match.arg(priorsizedistribution)
    prior <- dsizeprior(n=n,
		  type=priorsizedistribution,
		  mean.prior.size=mean.prior.size,
		  sd.prior.size=sd.prior.size,
		  mode.prior.sample.proportion=mode.prior.sample.proportion,
		  median.prior.sample.proportion=median.prior.sample.proportion,
		  median.prior.size=median.prior.size,
		  mode.prior.size=mode.prior.size,
		  quartiles.prior.size=quartiles.prior.size,
                  effective.prior.df=effective.prior.df,
                  maxN=maxN,
                  log=TRUE,
                  verbose=verbose)
    Cret <- .C("gplndisease",
              pop=as.integer(c(s,rep(0,prior$maxN-n))),
              dis=as.integer(c(dis,rep(0,prior$maxN-n))),
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
              N=as.integer(prior$N),
              maxN=as.integer(prior$maxN),
              sample=double(samplesize*dimsample),
              p0pos=double(K), p1pos=double(K),
              ppos=double(K), 
              lpriorm=as.double(prior$lprior),
              burnintheta=as.integer(burnintheta),
              verbose=as.integer(verbose), PACKAGE="size")
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
    #
    Cret$predictive.degree.count0<-Cret$nk0 / samplesize
    Cret$predictive.degree.count1<-Cret$nk1 / samplesize
    Cret$nk0 <- NULL
    Cret$nk1 <- NULL
    Cret$predictive.degree0<-Cret$p0pos
    Cret$p0pos<-NULL
    Cret$predictive.degree1<-Cret$p1pos
    Cret$p1pos<-NULL
    Cret$predictive.degree<-Cret$ppos
    Cret$ppos<-NULL
    #
  endrun <- burnin+interval*(samplesize-1)
  attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
  attr(Cret$sample, "class") <- "mcmc"
   mapfn <- function(x,lbound=min(x),ubound=max(x)){
     posdensN <- density(x, from=lbound, to=ubound)
     posdensN$x[which.max(posdensN$y)]
   }
   Cret$MAP <- apply(Cret$sample,2,mapfn)
   Cret$MAP["N"] <- mapfn(Cret$sample[,"N"],lbound=n,ubound=prior$maxN)
   Cret$MAP["disease"] <- mapfn(Cret$sample[,"disease"],lbound=0,ubound=1)
   Cret$maxN <- prior$maxN
   Cret$mode.prior.size <- prior$mode.prior.size
   Cret$effective.prior.df <- prior$effective.prior.df
   Cret$median.prior.size <- prior$median.prior.size
   Cret$mode.prior.sample.proportion <- prior$mode.prior.sample.proportion
   Cret$N <- prior$N
   Cret
}
