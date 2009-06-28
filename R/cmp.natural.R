cmp.natural <- function (mx,sdx,K=1000) 
{
#   natural is log(lambda), -nu
    options(warn = -1)
    result = optim(par=c(log(mx), 0), fn=function(p,mx,sdx,K) {
        j <- 0:K
        a <- dcmp(x=j, lambda=exp(p[1]), nu=exp(p[2]), err=0.000001)
	sqrt(abs(sum(a*j)^2-mx*mx)+abs(sum(a*j*j)-sdx*sdx-mx*mx))
    }, mx=mx, sdx=sdx, K=K, method = "Nelder-Mead",
       control=list(maxit=20000))
    options(warn = 0)
    lambda = exp(result$par[1])
    nu = exp(result$par[2])
    fit = c(lambda = lambda, nu = nu, result)
    return(fit)
}
cmp.mu <- function (p,K=1000) 
{
        j <- 0:K
        a <- dcmp(x=j, lambda=p[1], nu=p[2], err=0.000001)
	mu <- sum(a*j)
	sd <- sqrt(sum(a*j*j)-mu*mu)
    return(c(mu,sd))
}
