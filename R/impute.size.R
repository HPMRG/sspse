impute.size <- function(y, topcode = F, fix100k = F, na.rm = T, 
                   zero.to.na = T,
                   irule=c(0.25,        # 10 +/- 2.5
                           0.5,        # 20 +/- 2.5
                           0.5,        # 30 +/- 2.5
                           1.0,        # 40 +/- 2.5
                           2.5 ,        # 50 +/-  25
                           1.0,        # 60 +/- 2.5
                           1.0,        # 70 +/- 2.5
                           1.0,        # 80 +/- 2.5
                           1.0,        # 90 +/- 2.5
                           0,           # -  null coarsed +/- 0
                           10,           # 25 +/- 5
                           25,           # 75 +/- 5
                           1.25,        # 2500 +/- 125
                           1.25,        # 7500 +/- 125
                           0.1,         # 10000 +/- 1000
                           1)           # 15000 +/- 1000
                  )
{
#
# irule is the list of widths
#          1:4,6:9   the last digit of frag
#                5   y is not a multiple of 5
#               10   y is not a multiple of 10 (0)
#               11   y is ends in 25 (25)
#               12   y is ends in 75 (25)
#
#       irule=c(0.25,   # 10 +/- 2.5
#               0.25,   # 20 +/- 2.5
#               0.25,   # 30 +/- 2.5
#               0.25,   # 40 +/- 2.5
#               2.5 ,   # 50 +/-  25
#               0.25,   # 60 +/- 2.5
#               0.25,   # 70 +/- 2.5
#               0.25,   # 80 +/- 2.5
#               0.25,   # 90 +/- 2.5
#               0,      # -  null coarsed +/- 0
#               12.5,   # 25 +/- 12.5
#               12.5,   # 75 +/- 12.5
#               12.5,   # 2500 +/- 1250
#               12.5,   # 7500 +/- 1250
#               0.1,    # 10000 +/- 1000
#               1.0     # 15000 +/- 1000
# ASR irule
#
#       irule=c(0.1,    # 1
#               0.25,   # 1
#               0.25,   # 1
#               0.25,   # 1
#               0.25,   # 1
#               0.25,   # 1
#               0.25,   # 1
#               0.25,   # 1
#               0.25,   # 1
#               0,      # 1
#               2.5,    # 1
#               2.5,    # 1
#               12.5,   # 2500 +/- 1250
#               12.5,   # 7500 +/- 1250
#               0.1,    # 10000 +/- 1000
#               1.0     # 15000 +/- 1000
#
   yorig <- y      
#
#  Mark 0's in original
#
   if(zero.to.na) {
     y[y <= 0] <- NA
   }
#
#  Mark NAs in original
#
   if(na.rm) {
     yna <- is.na(y)
     y <- y[!yna]
   }
#
#  yrank <- floor(rank(y))
   yrank <- sort.list(sort.list(y))
#
   if(fix100k) {
#
#  first fix the 99999 people
#
    yup <- y[y > median(y)]
    p2 <- mean(yup < 99999)
    p1 <- 0.15
    xxx <- quantile(yup, c(p1, p2))
    a <- log((1 - p1)/(1 - p2))/log(xxx[2]/xxx[1])
    k <- xxx[1] * (1 - p1)^(1/a)
#
#   now sample
#
    nup <- y == 99999
    y[nup] <- round(k * (1 - (p2 + runif(sum(nup)) * (1 - p2)))^(-1/a))
   }
   if(!is.logical(topcode) | topcode) {
#
#  next fix the topcoding
#
    if(is.logical(topcode)) {
      topcode <- max(y)
    }
    yup <- y[y > median(y)]
    p2 <- mean(yup < topcode)
#   p2 <- 0.80
    p1 <- 0.15
    xxx <- quantile(yup, c(p1, p2))
    a <- log((1 - p1)/(1 - p2))/log(xxx[2]/xxx[1])
    k <- xxx[1] * (1 - p1)^(1/a)
#
#  now sample
#
    nup <- y >= topcode - 1e-09
    y[nup] <- round(k * (1 - (p2 + runif(sum(nup)) * (1 - p2)))^(-1/a))
#
#   Fix VERY non-pareto distributed values
#
    y[nup][y[nup] < topcode] <- topcode * 1.45
   }
#
#  create a table of the incomes
   tincome <- table(y)
#
#  y is the values and counts is the counts
#
   y <- as.numeric(names(tincome))
   counts <- as.numeric(tincome)
#
#  here
#
#  power10 is the baseline power of 10
#   power10 = 0 if y is not a power of 10
#   power10 = i if y is frag * 10^i
#  frag is the integer-multiple left over
#   y = frag * 10^power10
#
   power10 <- rep(0, length(y))
   for(i in (1:7)) {
    frag <- 10^(i) * floor(y/10^i)
    power10[y == frag] <- i
   }
   basis <- power10
   frag <- y/10^power10
#
#  r25 is the last digit of frag
#  basis is
#   1-9  the last digit of frag
#   10   y is not a multiple of 10
#   11   y is ends in 25
#
   r25 <- frag - 10 * floor(frag/10)
   basis[r25 == 5] <- 5
   r25 <- frag - 100 * floor(frag/100)
   basis[r25 == 25] <- 11
   basis[r25 == 75] <- 12
#
#  Pick up 2500, 7500
#
   basis[r25 == 25 & power10 > 1] <- 13
   basis[r25 == 75 & power10 > 1] <- 14
#
   r25 <- frag - 10 * floor(frag/10)
   for(i in c(1:4, 6:9)) {
    basis[r25 == i & basis < 11 & basis != 5] <- i
    basis[r25 == i & basis < 11 & power10 > 3] <- 15
   }
   basis[r25 == 5 & basis < 11 & power10 == 3] <- 16
#
#  Treat 50000 separately
#
   power10[r25 == 5 & power10 == 4] <- 3
#
#  Do not impute real 0's
#
   basis[y == 0] <- 10
#
#  basis[r25 == 5] <- 12
   basis[power10 == 0 & basis < 10 & basis != 5] <- 10
   imputew <- 10^(power10) * irule[basis]
#
#  xbasis is the list of incomes that do not need to be imputed
#
#  imputew is the width of the interval in dollars to impute on
#
   if(sum(basis == 10) > 0) {
    xbasis <- matrix(cbind(y, counts)[basis == 10,  ], ncol = 2)
    xbasis <- rep(xbasis[, 1], xbasis[, 2])
   }else{
    xbasis <- NULL
   }
#
#  isample is the table of values that do need to be imputed
#
#  The columns are: value, count, +/- width
#
#  imean is the set of values that do need to be imputed
#
   isample <- cbind(y, counts, imputew)[basis != 10,  ]
   imean <- rep(isample[, 1], isample[, 2])
#
#  iimput is the deviate to add to the value
#
   iimpute <- rep(isample[, 3], isample[, 2])*2*(runif(sum(isample[, 2])) - 0.5)
#  iimpute <- rep(isample[, 3], isample[, 2]) *
#            2 * (rep(1,sum(isample[, 2]))-0.5)
#
#  iimpute is the deviation to be added to each value
#  this wmi is the order og the values
#  the latter one is the imputed in their original order
   wmi <- order(c(xbasis, imean))
#  wmi <- c(xbasis, imean + round(iimpute))[wmi]
   wmi <- c(xbasis, imean + (iimpute))[wmi]
#
#  finally report them in the order they were entered
#  aaa <- order(yorig)
#  aaa <- cbind(yorig,wmi[yrank])[aaa,]
#  wmi[yrank]
#  wmi
   y <- yorig
   y[!yna] <- wmi[yrank]
   y
}
