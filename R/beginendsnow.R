beginsnow<-function(parallel=1, type="PVM", verbose=TRUE){
    ### snow is wrapper for MPI or PVM (mosix only has PVM)
    require(snow)
#
#   Start PVM if necessary
#
#   setDefaultClusterOptions(type="PVM")
    if(snow::getClusterOption("type")=="PVM"){
     if(verbose){
      cat("Engaging warp drive using PVM ...\n")
     }
     require(rpvm)
     PVM.running <- try(rpvm::.PVM.config(), silent = TRUE)
     if(inherits(PVM.running, "try-error")){
      hostfile <- paste(Sys.getenv("HOME"), "/.xpvm_hosts", sep = "")
      if(file.exists(hostfile)){
       rpvm::.PVM.start.pvmd(hostfile)
      }else{
       rpvm::.PVM.start.pvmd()
      }
      cat("no problem... PVM started by size...\n")
     }
    }else{
     if(verbose){
      cat("Engaging warp drive using MPI ...\n")
     }
    }
#
#   Start Cluster
#
    ### Snow commands to set up cluster
    cl <- snow::makeCluster(parallel)
    ### initialize parallel random number streams
    snow::clusterSetupRNG(cl)
    ### start each virtual machine with size library loaded
    snow::clusterEvalQ(cl, library(size))
#
#   Run the jobs with rpvm or Rmpi
#
    ### make sure that R has printed out console messages before go parallel
    flush.console()
    return(cl)
}
endsnow<-function(cl, verbose=TRUE){
    ### stop cluster and PVM (in case PVM is flakey)
    snow::stopCluster(cl)
    if(snow::getClusterOption("type")=="PVM") rpvm::.PVM.exit()
    invisible()
}
