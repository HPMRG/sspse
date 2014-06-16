cmp.natural <- function (mu,sig,K=1000,guess=NULL) 
{
# converts mean to natural
# natural is log(lambda), -nu
# natural is lambda, nu
    options(warn = -1)
    if(is.null(guess)){
      guess = c(2*log(mu+0.25), log(2))
    }
    result = optim(par=guess, fn=function(p,mu,sig,K) {
        j <- 0:K
        a <- dcmp(x=j, lambda=exp(p[1]), nu=exp(p[2]), err=0.000000001)
	sqrt(abs(sum(a*j)^2-mu*mu)+abs(sum(a*j*j)-sig*sig-mu*mu))
    }, mu=mu, sig=sig, K=K, method = "Nelder-Mead",
       control=list(maxit=20000))
    options(warn = 0)
    lambda = exp(result$par[1])
    nu = exp(result$par[2])
    fit = c(lambda = lambda, nu = nu, result)
    return(fit)
}
cmp.mu <- function (p,K=1000) 
{
# converts natural to mean
        j <- 0:K
        a <- dcmp(x=j, lambda=p[1], nu=p[2], err=0.000000001)
	mu <- sum(a*j)
	sd <- sqrt(sum(a*j*j)-mu*mu)
    return(c(mu,sd))
}
