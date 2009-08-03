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
