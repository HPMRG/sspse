beginparallel<-function(parallel=1, type="PVM", seed=NULL, verbose=TRUE){
    ### parallel is wrapper for MPI or PVM (mosix only has PVM)
#   require(snow)
    require(parallel)
#
#   Start PVM if necessary
#
#   setDefaultClusterOptions(type="PVM")
#   setDefaultClusterOptions(type="MPI")
#   if(snow::getClusterOption("type")=="PVM"){
    if(type=="PVM"){
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
    cl <- makeCluster(parallel,type=type)
    ### initialize parallel random number streams
    if(is.null(seed)){
     clusterSetRNGStream(cl)
    }else{
     clusterSetRNGStream(cl,iseed=seed)
    }
    ### start each virtual machine with size library loaded
    clusterEvalQ(cl, library(size))
#
#   Run the jobs with rpvm or Rmpi
#
    ### make sure that R has printed out console messages before go parallel
    flush.console()
    return(cl)
}
endparallel<-function(cl, type="PVM", verbose=TRUE){
    ### stop cluster and PVM (in case PVM is flakey)
    stopCluster(cl)
#   if(snow::getClusterOption("type")=="PVM") rpvm::.PVM.exit()
    if(type=="PVM") rpvm::.PVM.exit()
    invisible()
}
