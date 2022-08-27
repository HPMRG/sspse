#' @keywords internal
############### R-function: trapint ###############

# For performing trapezoidal integration.

# Last changed: 13 NOV 2003

.trapint <- function(xgrid,fgrid)
{
   ng <- length(xgrid)

   xvec <- xgrid[2:ng] - xgrid[1:(ng-1)]
   fvec <- fgrid[1:(ng-1)] + fgrid[2:ng]

   integ <- sum(xvec*fvec)/2

   return(integ)
}

############### End trapint ###############
