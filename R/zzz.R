######################################################################
# copyright (c) 2009, Krista J. Gile, University of Massachusetts - Amherst
#                     Mark S. Handcock, University of California - Los Angeles
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

.onLoad <-function(libname, pkgname){
	deducerNotLoaded <- try(.deducer == .jnull(),silent=TRUE)
	if(inherits(deducerNotLoaded,"try-error") || deducerNotLoaded)
		return(NULL)
	
  RDSAnalyst <- J("RDSAnalyst.RDSAnalyst")

  library.dynam("size", package=pkgname, lib.loc=libname)
  
  
  deducer.addMenuItem("Prior Distribution",,".getDialog('Prior Distribution')$run()","RDS Population")
  deducer.addMenuItem("Posterior Distribution",,".getDialog('Posterior Distribution')$run()","RDS Population")
# commented out marginal dialogs for now (JC 9/13/13)
#  deducer.addMenuItem("Marginal Posterior Distribution",,".getDialog('Marginal Posterior Distribution')$run()","RDS Population")
#  deducer.addMenuItem("Marginal Posterior Size Dist.",,".getDialog('Marginal Posterior Size Dist.')$run()","RDS Population")
  
  
  if(.windowsGUI){
	  winMenuAddItem("RDS Data",'Prior Distribution',"deducer('Prior Distribution')")
	  winMenuAddItem("RDS Data",'Posterior Distribution',"deducer('Posterior Distribution')")  
#	  winMenuAddItem("RDS Data",'Marginal Posterior Distribution',"deducer('Marginal Posterior Distribution')")  
#	  winMenuAddItem("RDS Data",'Marginal Posterior Size Dist.',"deducer('Marginal Posterior Size Dist.')")  
	  
  }
  
  else if(.jgr){	
	  if(RDSAnalyst$isPro()){
		  jgr.addMenuSeparator("RDS Population")
		  jgr.addMenuItem("RDS Population","Prior Distribution","deducer('Prior Distribution')")
		  jgr.addMenuItem("RDS Population","Posterior Distribution","deducer('Posterior Distribution')")
#		  jgr.addMenuItem("RDS Population","Marginal Posterior Distribution","deducer('Marginal Posterior Distribution')")
#		  jgr.addMenuItem("RDS Population","Marginal Posterior Size Dist.","deducer('Marginal Posterior Size Dist.')")
		  
	  }}
	
  .registerDialog("Prior Distribution", .makePriorDistribution)
  .registerDialog("Posterior Distribution", .makePosteriorDistribution)
#  .registerDialog("Marginal Posterior Distribution", .makeMarginalPosteriorDistribution)
#  .registerDialog("Marginal Posterior Size Dist.", .makeMarginalPosteriorSize)
  
}

.onAttach <- function(libname, pkgname){
  temp<-packageDescription("size")
  msg<-paste(temp$Package,": ",temp$Title,"\n",
      "Version ",temp$Version,
      " created on ",
      temp$Date,".\n", sep="")
  msg<-paste(msg,"copyright (c) 2009, Krista J. Gile, University of Massachusetts - Amherst\n",
"                    Mark S. Handcock, University of California - Los Angeles\n",sep="")
  msg<-paste(msg,'For citation information, type citation("size").\n')
  msg<-paste(msg,'Type help("size-package") to get started.\n')
  packageStartupMessage(msg)
}
