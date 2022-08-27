poscmpwp<-function(s,s2=NULL,rc=rep(FALSE,length=length(s2)),maxN=NULL,
                  K=2*max(c(s,s2)), n=length(s), n2=length(s2),
                  mean.prior.visibility=7, sd.prior.visibility=3,
                  df.mean.prior.visibility=1, df.sd.prior.visibility=5,
                  beta0.mean.prior=-3, beta1.mean.prior=0,
                  beta0.sd.prior=10, beta1.sd.prior=10,
                  mem.optimism.prior=1, df.mem.optimism.prior=5, 
                  mem.scale.prior=1, df.mem.scale.prior=5, 
		  mem.overdispersion=5,
                  muproposal=0.1, 
                  nuproposal=0.15, 
                  beta0proposal=0.1, beta1proposal=0.001,
                  memmuproposal=0.1, memscaleproposal=0.15,
                  visibility=TRUE,
                  Np=0,
                  samplesize=10,burnin=0,interval=1,burnintheta=500,burninbeta=20,
                  priorsizedistribution=c("beta","nbinom","pln","flat","supplied"),
                  mean.prior.size=NULL, sd.prior.size=NULL,
                  mode.prior.sample.proportion=NULL,
                  median.prior.sample.proportion=NULL,
                  median.prior.size=NULL,
                  mode.prior.size=NULL,
                  quartiles.prior.size=NULL,
                  effective.prior.df=1,
                  alpha=NULL,
                  seed=NULL,
                  maxbeta=120,
                  supplied=list(maxN=maxN),
              num.recruits=NULL,
              recruit.times=NULL,
              num.recruits2=NULL,
              recruit.times2=NULL,
              max.coupons=NULL,
                  verbose=TRUE){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    if(!is.null(seed))  set.seed(as.integer(seed))
    #
    # Cap the maximum visibility to K
    #
    # s[s>K] <- K
    #
    # Transform observed mean parametrization to log-normal
    # parametrization
    #
    lnlam <- log(mean.prior.visibility)
    nu <- sd.prior.visibility*sd.prior.visibility/mean.prior.visibility
    #
    if(visibility){
      dimsample <- 5+Np+4
    }else{
      dimsample <- 5+Np
    }
    #
    # Determine if we are in the two-sample case or the one-sample
    n1 = n
    if(!is.null(s2)){
     n0 = sum(rc)
     n = n1 + n2 - n0 # The number of unique people seen
    }
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
                  maxbeta=maxbeta,
                  log=TRUE,
                  supplied=supplied,
                  verbose=verbose)
    message(sprintf("Maximum population size set to %d.\n",prior$maxN),appendLF=FALSE)
    if(!is.null(s2)){
     if(visibility){
      cat(sprintf("Using Capture-recapture Weighted Negative Binomial measurement error model with K = %d.\n",K))
      Cret <- .C("gcmpwpvis2",
              pop12=as.integer(c(s, s2[!rc], rep(0,prior$maxN-length(s)-length(s2[!rc])))),
              pop21=as.integer(c(s2,sum(s)-sum(s2[rc]), rep(0,prior$maxN-length(s2)-1))),
              K=as.integer(K),
              n1=as.integer(n1),
              n2=as.integer(n2),
              n0=as.integer(n0),
              samplesize=as.integer(samplesize),
              burnin=as.integer(burnin),
              interval=as.integer(interval),
              mu=as.double(mean.prior.visibility), df.mean.prior.visibility=as.double(df.mean.prior.visibility),
              sigma=as.double(sd.prior.visibility), df.sd.prior.visibility=as.double(df.sd.prior.visibility),
              lnlam=as.double(lnlam), nu=as.double(nu),
              beta0.mean.prior=as.double(beta0.mean.prior), beta0.sd.prior=as.double(beta0.sd.prior),
              beta1.mean.prior=as.double(beta1.mean.prior), beta1.sd.prior=as.double(beta1.sd.prior),
              mem.optimism.prior=as.double(log(mem.optimism.prior)), df.mem.optimism.prior=as.double(df.mem.optimism.prior),
              mem.scale.prior=as.double(mem.scale.prior^2), df.mem.scale.prior=as.double(df.mem.scale.prior),
              mem.overdispersion=as.double(mem.overdispersion),
              Np=as.integer(Np),
              srd=as.integer(s),
              numrec=as.integer(num.recruits),
              rectime=as.double(recruit.times),
              srd2=as.integer(s2),
              numrec2=as.integer(num.recruits2),
              rectime2=as.double(recruit.times2),
              rc=as.integer(rc),
              maxcoupons=as.integer(max.coupons),
              muproposal=as.double(muproposal),
              nuproposal=as.double(nuproposal),
              beta0proposal=as.double(beta0proposal), beta1proposal=as.double(beta1proposal),
              memmuproposal=as.double(memmuproposal), memscaleproposal=as.double(memscaleproposal),
              N=as.integer(prior$N),
              maxN=as.integer(prior$maxN),
              sample=double(samplesize*dimsample),
              vsample=integer(samplesize*n1),
              vsample2=integer(samplesize*n2),
              posu=double(K),
              posd=double(10*K),
              lpriorm=as.double(prior$lprior),
              burnintheta=as.integer(burnintheta),
              burninbeta=as.integer(burninbeta),
              verbose=as.integer(verbose), PACKAGE="sspse")
     }else{
      cat(sprintf("Using Capture-recapture non-measurement error model with K = %d.\n",K))
       s[ s>K] <- K
      s2[s2>K] <- K
       s[ s<1] <- 1
      s2[s2<1] <- 1
      nk=tabulate(c(s,s2[!rc]),nbins=K)
      Cret <- .C("gcmp2",
              pop12=as.integer(c(s, s2[!rc], rep(0,prior$maxN-length(s)-length(s2[!rc])))),
              pop21=as.integer(c(s2,sum(s)-sum(s2[rc]), rep(0,prior$maxN-length(s2)-1))),
              nk=as.integer(nk),
              K=as.integer(K),
              n1=as.integer(n1),
              n2=as.integer(n2),
              n0=as.integer(n0),
              samplesize=as.integer(samplesize),
              burnin=as.integer(burnin),
              interval=as.integer(interval),
              mu=as.double(mean.prior.visibility), df.mean.prior.visibility=as.double(df.mean.prior.visibility),
              sigma=as.double(sd.prior.visibility), df.sd.prior.visibility=as.double(df.sd.prior.visibility),
              lnlam=as.double(lnlam), nu=as.double(nu),
              Np=as.integer(Np),
              muproposal=as.double(muproposal),
              nuproposal=as.double(nuproposal),
              N=as.integer(prior$N),
              maxN=as.integer(prior$maxN),
              sample=double(samplesize*dimsample),
              posu=double(K),
              lpriorm=as.double(prior$lprior),
              burnintheta=as.integer(burnintheta),
              verbose=as.integer(verbose), PACKAGE="sspse")
     }
    }else{
     if(visibility){
      cat(sprintf("Using Weighted Negative Binomial measurement error model with K = %d.\n",K))
      Cret <- .C("gcmpwpvis",
              pop=as.integer(c(s,rep(0,prior$maxN-n1))),
              K=as.integer(K),
              n=as.integer(n1),
              samplesize=as.integer(samplesize),
              burnin=as.integer(burnin),
              interval=as.integer(interval),
              mu=as.double(mean.prior.visibility), df.mean.prior.visibility=as.double(df.mean.prior.visibility),
              sigma=as.double(sd.prior.visibility), df.sd.prior.visibility=as.double(df.sd.prior.visibility),
              lnlam=as.double(lnlam), nu=as.double(nu),
              beta0.mean.prior=as.double(beta0.mean.prior), beta0.sd.prior=as.double(beta0.sd.prior),
              beta1.mean.prior=as.double(beta1.mean.prior), beta1.sd.prior=as.double(beta1.sd.prior),
              mem.optimism.prior=as.double(log(mem.optimism.prior)), df.mem.optimism.prior=as.double(df.mem.optimism.prior),
              mem.scale.prior=as.double(mem.scale.prior^2), df.mem.scale.prior=as.double(df.mem.scale.prior),
              mem.overdispersion=as.double(mem.overdispersion),
              Np=as.integer(Np),
              srd=as.integer(s),
              numrec=as.integer(num.recruits),
              rectime=as.double(recruit.times),
              maxcoupons=as.integer(max.coupons),
              muproposal=as.double(muproposal),
              nuproposal=as.double(nuproposal),
              beta0proposal=as.double(beta0proposal), beta1proposal=as.double(beta1proposal),
              memmuproposal=as.double(memmuproposal), memscaleproposal=as.double(memscaleproposal),
              N=as.integer(prior$N),
              maxN=as.integer(prior$maxN),
              sample=double(samplesize*dimsample),
              vsample=integer(samplesize*n1),
              posu=double(K),
              posd=double(10*K),
              lpriorm=as.double(prior$lprior),
              burnintheta=as.integer(burnintheta),
              burninbeta=as.integer(burninbeta),
              verbose=as.integer(verbose), PACKAGE="sspse")
     }else{
      cat(sprintf("Using non-measurement error model with K = %d.\n",K))
      s[s>K] <- K
      s[s<1] <- 1
      Cret <- .C("gcmp",
              pop=as.integer(c(s,rep(0,prior$maxN-n1))),
              K=as.integer(K),
              n=as.integer(n1),
              samplesize=as.integer(samplesize),
              burnin=as.integer(burnin),
              interval=as.integer(interval),
              mu=as.double(mean.prior.visibility), df.mean.prior.visibility=as.double(df.mean.prior.visibility),
              sigma=as.double(sd.prior.visibility), df.sd.prior.visibility=as.double(df.sd.prior.visibility),
              lnlam=as.double(lnlam), nu=as.double(nu),
              Np=as.integer(Np),
              muproposal=as.double(muproposal),
              nuproposal=as.double(nuproposal),
              N=as.integer(prior$N),
              maxN=as.integer(prior$maxN),
              sample=double(samplesize*dimsample),
              posu=double(K),
              lpriorm=as.double(prior$lprior),
              burnintheta=as.integer(burnintheta),
              verbose=as.integer(verbose), PACKAGE="sspse")
     }
    }
    Cret$sample<-matrix(Cret$sample,nrow=samplesize,ncol=dimsample,byrow=TRUE)
    degnames <- NULL
    if(Np>0){degnames <- c(degnames,paste("pdeg",1:Np,sep=""))}
    if(visibility){
     colnamessample <- c("N","mu","sigma","visibility1","totalsize","beta0","beta1","mem.optimism","mem.scale")
     Cret$vsample<-matrix(Cret$vsample,nrow=samplesize,ncol=n1,byrow=TRUE)
     colnames(Cret$vsample) <- 1:n1
     if(!is.null(s2)){
      Cret$vsample2<-matrix(Cret$vsample2,nrow=samplesize,ncol=n2,byrow=TRUE)
      colnames(Cret$vsample2) <- 1:n2
     }
     max.mu <- 3*median(Cret$vsample)
     if(length(degnames)>0){
      colnames(Cret$sample) <- c(colnamessample, degnames)
     }else{
      colnames(Cret$sample) <- colnamessample
     }
     if(stats::var(Cret$sample[,"mem.optimism"])>1e-8){
       Cret$sample[,"mem.optimism"] <- exp(Cret$sample[,"mem.optimism"])
     }else{
       Cret$sample <- Cret$sample[,-match("mem.optimism",colnames(Cret$sample))]
     }
     #
     # Transform WP parametrization to mean
     # parametrization (for the mean visibility value)
     #
     Cret$sample[,"mem.scale"] <- sqrt(Cret$sample[,"mem.scale"])
#    a <- Cret$sample[,c("mem.optimism","mem.scale")]
#    a[,2] <- sqrt(a[,1]*mean(Cret$vsample)*a[,2])
#    colnames(a) <- c("mem.optimism","mem.sigma")
#    nas <- apply(a,1,function(x){any(is.na(x))})
#    if(!all(nas)){
#     inas <- sample(seq_along(nas)[!nas],size=sum(nas),replace=TRUE)
#     a[nas,] <- a[inas,]
##    Cret$sample[,c("mem.optimism","mem.scale")] <- a
#     Cret$sample <- cbind(Cret$sample,a[,2])
#     colnames(Cret$sample)[ncol(Cret$sample)] <- "mem.sigma"
#    }
    }else{
     colnamessample <- c("N","mu","sigma","visibility1","totalsize")
     max.mu <- 2*mean.prior.visibility
     if(length(degnames)>0){
      colnames(Cret$sample) <- c(colnamessample, degnames)
     }else{
      colnames(Cret$sample) <- colnamessample
     }
    }
   if(F){
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
    a <- t(apply(Cret$sample[,c("mu","sigma")],1,cmp.to.mu.sd, max.mu=max.mu))
    nas <- apply(a,1,function(x){any(is.na(x))})
    if(!all(nas)){
     inas <- sample(seq_along(nas)[!nas],size=sum(nas),replace=TRUE)
     a[nas,] <- a[inas,]
#    Cret$sample[,c("mu","sigma")] <- t(apply(Cret$sample[,c("mu","sigma")],1,cmp.to.mu.sd,max.mu=5*mean.prior.visibility)))
     Cret$sample[,c("mu","sigma")] <- a
    }else{
      warning(paste("All the lambda and nu parameters are extreme. The mean and sigma are on the natural scale."), call. = FALSE)
    }
   }
#   # PLN mean
#   Cret$sample <- cbind(Cret$sample,Cret$sample[,c("mem.optimism")])
#   colnames(Cret$sample)[ncol(Cret$sample)] <- c("mem.visibility.mean")
#   mean.visibility <- sum(seq(along=Cret$nk)*Cret$posu)
#   print(mean.visibility)
#   Cret$sample[,"mem.visibility.mean"] <- exp(log(mean.visibility)+Cret$sample[,"mem.optimism"]+0.5*Cret$sample[,"mem.scale"])
    #
#   Cret$Nk<-Cret$nk/sum(Cret$nk)
    Cret$predictive.visibility.count<-Cret$nk / samplesize
    Cret$nk<-NULL
    Cret$predictive.visibility<-Cret$posu
    Cret$predictive.degree<-Cret$posd
    Cret$posu<-NULL
    endrun <- burnin+interval*(samplesize-1)
    attr(Cret$sample, "mcpar") <- c(burnin+1, endrun, interval)
    attr(Cret$sample, "class") <- "mcmc"
    ### compute modes of posterior samples,Maximum A Posterior (MAP) values N, mu, sigma, visibility1
    Cret$MAP <- apply(Cret$sample,2,mode.density)
    Cret$MAP["N"] <- mode.density(Cret$sample[,"N"],lbound=n,ubound=prior$maxN)
    if(!is.null(s2)){Cret$n <- Cret$n1 +  Cret$n2 - Cret$n0}
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
