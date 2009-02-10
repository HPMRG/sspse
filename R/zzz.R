######################################################################
# copyright (c) 2009, Krista J. Gile, Nuffield College, Oxford,
#                     Mark S. Handcock, University of Washington
# 
# For license and citation information see
#    http://statnet.org/attribution
#
# We have invested a lot of time and effort in creating 'statnet',
# for use by other researchers. We require that the attributions
# in the software are retained (even if only pieces of it are used),
# and that there is attribution when the package is loaded (e.g., via
# "library" or "require"). 
######################################################################
# File name: zzz.R
######################################################################
#
# .First.lib is run when the package is loaded.
#
######################################################################

.First.lib <- function(lib, pkg){
  library.dynam("size", pkg, lib)
  DESCpath <- file.path(system.file(package="size"), "DESCRIPTION")
  info <- read.dcf(DESCpath)
  cat('\nsize:', info[,"Title"], 
      '\nVersion', info[,"Version"], 'created on', info[,"Date"], '\n') 
  
  cat(paste("copyright (c) 2009, Krista J. Gile, Nuffield College, Oxford\n",
"                    Mark S. Handcock, University of Washington\n",sep=""))
  cat('Type help(package="size") to get started.\n\n')
  cat('Based on "statnet" project software (http://statnet.org).\n',
      'For license and citation information see http://statnet.org/attribution\n',
      'or type citation("size").\n')
}

.Last.lib <- function(libpath){
  library.dynam.unload("size",libpath)
}
