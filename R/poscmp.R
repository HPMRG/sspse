poscmp<-function(s,maxN=4*length(s),
                  K=2*max(s), nk=tabulate(s,nbin=K), n=length(s),
		  N=maxN/2,
                  mean.prior.degree=7, sd.prior.degree=3,
                  df.mean.prior=1, df.sd.prior=5,
                  muproposal=0.25, 
                  sigmaproposal=0.15, 
                  Np=0,
                  samplesize=10,burnin=0,interval=1,burnintheta=500,
		  mean.prior.size=N, sd.prior.size=N,
                  seed=NULL,
                  verbose=TRUE){
#		  mean.prior.size=N, sd.prior.size=ceiling(sqrt(mean.prior.size)*3),
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    out <- cmp.natural(mean.prior.degree, sd.prior.degree)
    mu0 <- log(out$lambda)
    sigma0 <- out$nu
#print(c(mu0,sigma0))
    if(!is.null(seed))  set.seed(as.integer(seed))
    dimsample <- 4+Np
    lpriorm <- dnbinommu(x=n+(1:maxN)-1,
                         mu=mean.prior.size, sd=sd.prior.size,
			 log=TRUE)
    Cret <- .C("gcmp",
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
    Cret$sample[,"mu"] <- exp(Cret$sample[,"mu"])
    Cret$sample <- cbind(Cret$sample,Cret$sample[,c("mu","sigma")])
    colnames(Cret$sample)[ncol(Cret$sample)-(1:0)] <- c("lambda","nu")
    # Transform to mean value parametrization 
    Cret$sample[,c("mu","sigma")] <- t(apply(Cret$sample[,c("mu","sigma")],1,cmp.mu))
    #
    Cret$Nk<-Cret$nk/sum(Cret$nk)
    endrun <- burnin+interval*(samplesize-1)
    attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    mapfn <- function(x,lbound=min(x),ubound=max(x)){
      posdensN <- density(x, from=lbound, to=ubound)
      posdensN$x[which.max(posdensN$y)]
    }
    ### compute modes of posterior samples,Maximum A Posterior (MAP) values N, mu, sigma, degree1
    Cret$MAP <- apply(Cret$sample,2,mapfn)
    Cret$MAP["N"] <- mapfn(Cret$sample[,"N"],lbound=n,ubound=maxN)
    Cret
}

rcmp.lambda <- function(n, lambda, nu, err=0.00001, K=100){
  # Perform argument checking
  if (lambda < 0 || nu < 0)
	stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
  if (n < 0 || n != floor(n))
	stop("Invalid number of draws");
  out <- .C("rcmp",
            x=integer(n),
            lambda=as.double(lambda),
            nu=as.double(nu),
            n=as.integer(n),
            K=as.integer(K),
            err=as.double(err),
            PACKAGE="size")
   return(out$x)
}
rcmp <- function(n, mu, sig, err=0.00001, K=100){
  out <- cmp.natural(mu,sig)
  # Perform argument checking
  if (out$lambda < 0 || out$nu < 0)
	stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
  if (n < 0 || n != floor(n))
	stop("Invalid number of draws");
  out <- .C("rcmp",
            x=integer(n),
            lambda=as.double(out$lambda),
            nu=as.double(out$nu),
            n=as.integer(n),
            K=as.integer(K),
            err=as.double(err),
            PACKAGE="size")
   return(out$x)
}
dcmp <- function(x, lambda, nu, err=0.00001, log=FALSE){
  # Perform argument checking
  if (lambda < 0 || nu < 0)
	stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
  if (any(x < 0 || x != floor(x)))
		return (0);
  out <- .C("dcmp",
            x=as.integer(x),
            lambda=as.double(lambda),
            nu=as.double(nu),
            n=as.integer(length(x)),
            err=as.double(err),
            give_log=as.integer(log),
            val=double(length(x)),
            PACKAGE="size")
   return(out$val)
}
dcmp.mu <- function(x, mu, sig, err=0.00001, log=FALSE){
  # Perform argument checking
  if (mu < 0 || sig < 0)
	stop("Invalid arguments, only defined for mu >= 0, sd >= 0");
  if (any(x < 0 || x != floor(x)))
		return (0);
  out <- cmp.natural(mu,sig)
  out <- .C("dcmp",
            x=as.integer(x),
            lambda=as.double(out$lambda),
            nu=as.double(out$nu),
            n=as.integer(length(x)),
            err=as.double(err),
            give_log=as.integer(log),
            val=double(length(x)),
            PACKAGE="size")
   return(out$val)
}
cmp.compute.z = function(lambda, nu, error = 0.001, log=FALSE)
{
    if (missing(lambda)) {
      stop("The mean of the Negative Binomial must be specified.")
    }
    if (missing(nu)) {
      stop("The standard deviation of the Negative Binomial must be specified.")
    }
#  return(out$val)
   out <- .C("vzcmp",
            lambda=as.double(lambda),
            nu=as.double(nu),
            err=as.double(error),
            give_log=as.integer(log),
            out=as.double(error),
            PACKAGE="size")
   return(out$out)
}
dcomsize = function(x, lambda, nu, z = NULL)
{
	# Perform argument checking
	if (lambda < 0 || nu < 0)
		stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
	if (x < 0 || x != floor(x))
		return (0);
	if (is.null(z) || z <= 0)
		z = cmp.compute.z(lambda, nu);
	
	# Return pmf
	return ((lambda ^ x) * ((factorial(x)) ^ -nu) / z);
}
