unposN<-function(mu,rho,N,K,s,n=tabulate(s,nbin=K),M=100000){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    if(length(s)>N){print("Error - observed unit outside range")}
    prob=rep(0,K)
    Cret <- .C("bnw_NC",
              N=as.integer(N),
              K=as.integer(K),
              n=as.integer(length(s)),
              s=as.integer(s),
              nk=as.integer(n),
              Nk=as.integer(prob),
              prob=as.double(prob),
              mu=as.double(mu),
              rho=as.double(rho),
              M=as.integer(M),
              unpos=as.double(0))
    Cret$unpos
    Cret
}
