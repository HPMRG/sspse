cmpmle <-function(xv,xp,cutoff=1,cutabove=1000,guess=c(7,3)){
 xp <- xp[xv >= cutoff]
 xv <- xv[xv >= cutoff]
 xp <- xp[xv <= cutabove]
 xv <- xv[xv <= cutabove]
 xp <- length(xv)*xp/sum(xp)
 guess <- cmp.to.natural(mu=guess[1],sigma=guess[2])
 guess <- c(log(guess$lambda), log(guess$nu))
  aaa <- stats::optim(par=guess,fn=allcmp,
   method="BFGS",
   hessian=FALSE,control=list(fnscale=-10),
   xv=xv,xp=xp,cutoff=cutoff,cutabove=cutabove)
  names(aaa$par) <- c("CMP lambda","CMP nu")
  exp(aaa$par)
}
nb0mle <-function(xv,xp,cutoff=1,cutabove=1000,guess=c(5,0.1)){
 xp <- xp[xv >= cutoff]
 xv <- xv[xv >= cutoff]
 xp <- xp[xv <= cutabove]
 xv <- xv[xv <= cutabove]
 xp <- length(xv)*xp/sum(xp)
  aaa <- stats::optim(par=guess,fn=allnb0,
   method="BFGS",
   hessian=FALSE,control=list(fnscale=-10),
   xv=xv,xp=xp,cutoff=cutoff,cutabove=cutabove)
  aaa <- nbmean(aaa)
  names(aaa) <- c("gamma mean","gamma s.d.")
  aaa
}
allnb0 <- function(v,xv,xp,cutoff=1,cutabove=1000){
 n <- sum(xp)
 if(v[1]<=0 | v[1]>=15000 | v[2]<=0 | v[2]>=1){
  out <- -10^10
 }else{
  if(cutabove<1000){
   cprob <- sum(stats::dnbinom(x=cutoff:cutabove, size=v[1]*v[2], prob=v[2]))
  }else{
   cprob <- stats::pnbinom(q=cutoff-1, size=v[1]*v[2], prob=v[2],lower.tail=FALSE)
  }
  lpgp <- stats::dnbinom(x=xv, size=v[1]*v[2], prob=v[2],log=TRUE)
  out <- sum(xp*lpgp)-n*log(cprob)
  if(is.infinite(out)){out <- -10^10}
  if(is.na(out)){out <- -10^10}
 }
 out
}
nbmean <- function(aaa=NULL,theta=NULL){
#
# Input is (expected number of successes, prob of success)
#           = alpha*(1+beta), 1/(1+beta)
#
# E(K) = E(success) + E(failure) = size + E(X) = alpha + alpha*beta
# V(K) = V(failure) = V(X) = alpha*beta*beta
#
# If theta ~ Gamma(alpha, beta)
# then X ~ NB(alpha, beta)
#
# alpha = size
# beta = scale = (1-prob)/prob
#
#  mean gamma = alpha*beta
#  var  gamma = alpha*beta*beta
#
#  mean degree = E(K) = alpha + alpha*beta
#  var  degree = V(K) = alpha*beta*beta
#
# Output is (gamma mean, gamma s.d.)
#
 if(!is.null(theta)){aaa$theta <- theta}
 if(is.null(aaa$theta)){aaa$theta <- aaa$par}
#if(aaa$theta[1]>1){aaa$theta[1] <- 1/aaa$theta[1]}
#
#gmu <- aaa$theta[1]*aaa$theta[2]
#var <- mu + mu*mu/aaa$theta[2]
#
 prob <- aaa$theta[2]
 alpha <- prob*aaa$theta[1]
 beta <- (1-prob)/prob
 mu <- alpha*beta
 var <- alpha*beta*beta
 bbb <- c(mu,sqrt(var))
 names(bbb) <- c("gamma mean","gamma s.d.")
 bbb
}
allcmp <- function(v,xv,xp,cutoff=1,cutabove=1000){
 n <- sum(xp)
 if(cutabove<1000){
  cprob <- sum(dcmp.natural(v=exp(v),x=(cutoff:cutabove)))
 }else{
  if(cutoff>0){
   cprob <- 1 - sum(dcmp.natural(v=exp(v),x=0:(cutoff-1)))
  }else{
   cprob <- 1
  }
 }
#
 out <- sum(xp*ldcmp.natural(v=exp(v),x=xv))-n*log(cprob)
 if(is.infinite(out)){out <- -10^10}
 if(is.na(out)){out <- -10^10}
 out
}
#
# Compute the CMP PMF
#
ldcmp.natural <- function(v,x,cutoff=0){
 log(dcmp.natural(v,x,cutoff=cutoff))
}
dcmp.natural <- function(v, x, cutoff=0, err=0.000000001, log=FALSE){
  # compute PMF from natural parameters
  # Perform argument checking
  if (v[1] < 0 | v[2] < 0)
	stop("Invalid arguments, only defined for mu >= 0, sd >= 0");
  y <- .C("dcmp",
            x=as.integer(x),
            lambda=as.double(v[1]),
            nu=as.double(v[2]),
            n=as.integer(length(x)),
            err=as.double(err),
            give_log=as.integer(log),
            val=double(length(x)),
            PACKAGE="sspse")$val
  if (cutoff > 0) {
    c0 <- 1 - sum(dcmp.natural(v = v, x = (0:(cutoff - 1)), cutoff=0))
    y <- y/c0
  }
  return(y)
}
dcmp <- function(x, lambda, nu, err=0.000000001, log=FALSE){
  # compute PMF from natural parameters (with different parsing arguments)
  # Perform argument checking
  if (lambda < 0 | nu < 0)
	stop("Invalid arguments, only defined for lambda >= 0, nu >= 0");
  if (any(x < 0 | (x != floor(x))))
		return (0);
  out <- .C("dcmp",
            x=as.integer(x),
            lambda=as.double(lambda),
            nu=as.double(nu),
            n=as.integer(length(x)),
            err=as.double(err),
            give_log=as.integer(log),
            val=double(length(x)),
            PACKAGE="sspse")
   return(out$val)
}
