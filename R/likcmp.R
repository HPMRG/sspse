likcmp<-function(s,maxN=NULL,
                  K=2*max(s), nk=tabulate(s,nbins=K), n=length(s),
                  mean.prior.visibility=7, sd.prior.visibility=3,
                  df.mean.prior.visibility=1, df.sd.prior.visibility=5,
                  muproposal=0.1, 
                  sigmaproposal=0.15, 
                  Np=0,
                  samplesize=10,warmup=0,interval=1,warmuptheta=500,
		  priorsizedistribution=c("flat","beta","nbinom","pln"),
		  mean.prior.size=NULL, sd.prior.size=NULL,
		  mode.prior.sample.proportion=0.5,
		  median.prior.sample.proportion=NULL,
		  median.prior.size=NULL,
		  mode.prior.size=NULL,
		  quartiles.prior.size=NULL,
		  effective.prior.df=1,
		  alpha=NULL,
                  seed=NULL,
                  verbose=TRUE){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    if(!is.null(seed))  set.seed(as.integer(seed))
    #
    # Cap the maximum degree to K
    #
    s[s>K] <- K
    #
    # Transform observed mean parametrization to log-normal
    # parametrization
    #
    out <- cmp.to.natural(mu=mean.prior.visibility, sigma=sd.prior.visibility)
    mu <- log(out$lambda)
    sigma <- out$nu
    #
    dimsample <- 4+Np
    #
    priorsizedistribution=match.arg(priorsizedistribution)
    prior <- dsizeprior(n=n,
		  type=priorsizedistribution,
		  sd.prior.size=sd.prior.size,
		  mode.prior.sample.proportion=mode.prior.sample.proportion,
		  median.prior.sample.proportion=median.prior.sample.proportion,
		  median.prior.size=median.prior.size,
		  mode.prior.size=mode.prior.size,
		  mean.prior.size=mean.prior.size,
		  quartiles.prior.size=quartiles.prior.size,
                  effective.prior.df=effective.prior.df,
                  alpha=alpha,
                  maxN=maxN,
                  log=TRUE,
                  verbose=verbose)
    Cret <- .C("lcmp",
              pop=as.integer(c(s,rep(0,prior$maxN-n))),
              nk=as.integer(nk),
              K=as.integer(K),
              n=as.integer(n),
              samplesize=as.integer(samplesize),
              warmup=as.integer(warmup),
              interval=as.integer(interval),
              mu=as.double(mu), df.mean.prior.visibility=as.double(df.mean.prior.visibility),
              sigma=as.double(sigma), df.sd.prior.visibility=as.double(df.sd.prior.visibility),
              Np=as.integer(Np),
              muproposal=as.double(muproposal),
              sigmaproposal=as.double(sigmaproposal),
              N=as.integer(prior$N),
              maxN=as.integer(prior$maxN),
              sample=double(samplesize*dimsample),
              ppos=double(K),
              lpriorm=as.double(prior$lprior),
              warmuptheta=as.integer(warmuptheta),
              verbose=as.integer(verbose), PACKAGE="sspse")
    Cret$sample<-matrix(Cret$sample,nrow=samplesize,ncol=dimsample,byrow=TRUE)
    degnames <- NULL
    if(Np>0){degnames <- c(degnames,paste("pdeg",1:Np,sep=""))}
    colnamessample <- c("N","mu","sigma","degree1")
    if(length(degnames)>0){
     colnames(Cret$sample) <- c(colnamessample, degnames)
    }else{
     colnames(Cret$sample) <- colnamessample
    }
    #
    # Transform observed mean parametrization to log-normal
    # parametrization
    #
    # Expectation and s.d. of normal from log-normal
    #
    Cret$sample[,"mu"] <- exp(Cret$sample[,"mu"])
    Cret$sample <- cbind(Cret$sample,Cret$sample[,c("mu","sigma")])
    colnames(Cret$sample)[ncol(Cret$sample)-(1:0)] <- c("lambda","nu")
    # Transform to mean value parametrization 
    Cret$sample[,c("mu","sigma")] <- t(apply(Cret$sample[,c("mu","sigma")],1,cmp.to.mu.sd,max.mu=5*mean.prior.visibility))
    #
#   Cret$Nk<-Cret$nk/sum(Cret$nk)
    Cret$predictive.visibility.count<-Cret$nk / samplesize
    Cret$nk<-NULL
    Cret$predictive.visibility<-Cret$ppos
    Cret$ppos<-NULL
    endrun <- warmup+interval*(samplesize-1)
    attr(Cret$sample, "mcpar") <- c(warmup+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    ### compute modes of posterior samples,Maximum A Posterior (MAP) values N, mu, sigma, degree1
    Cret$MAP <- apply(Cret$sample,2,mode.density)
    Cret$MAP["N"] <- mode.density(Cret$sample[,"N"],lbound=n,ubound=prior$maxN)
#
#   Cret$MSE <- c(((prior$x-mean.prior.visibility)^2)*prior$lprior/sum(prior$lprior),mean((Cret$sample[,"N"]-mean.prior.visibility)^2))
    Cret$maxN <- prior$maxN
    Cret$quartiles.prior.size <- prior$quartiles.prior.size
    Cret$mode.prior.size <- prior$mode.prior.size
    Cret$mean.prior.size <- prior$mean.prior.size
    Cret$effective.prior.df <- prior$effective.prior.df
    Cret$median.prior.size <- prior$median.prior.size
    Cret$mode.prior.sample.proportion <- prior$mode.prior.sample.proportion
    Cret$N <- prior$N
    Cret
}
