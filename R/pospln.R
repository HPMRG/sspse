pospln<-function(s,mu0=2,kappa0=2/0.06,sigma20=0.06,df0=5,
                 sigma2proposal=0.002, 
                 samplesize=1000,burnin=0,interval=1,
                 seed=NULL,
                 verbose=TRUE){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    if(!is.null(seed))  set.seed(as.integer(seed))
    musample <- rep(0,samplesize)
    sigmasample <- rep(0,samplesize)
    Cret <- .C("MetropolisHastings",
              s=as.integer(s),
              mu0=as.double(mu0), kappa0=as.double(kappa0),
              sigma20=as.double(sigma20), df0=as.double(df0),
              sigma2proposal=as.double(sigma2proposal),
              N=as.integer(length(s)),
              musample=as.double(musample),
              sigmasample=as.double(sigmasample),
              nsteps=as.integer(samplesize),
              staken=as.integer(0),
              burnin=as.integer(burnin),
              interval=as.integer(interval),
              fVerbose=as.integer(verbose))
    Cret$sample <- cbind(Cret$musample, Cret$sigmasample)
    colnames(Cret$sample) <- c("mu","sigma")
    endrun <- burnin+interval*(samplesize-1)
    attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    Cret
}
