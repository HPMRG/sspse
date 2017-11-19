#' @keywords internal
beginparallel<-function(parallel=1, type="PSOCK", seed=NULL, packagenames=c("sspse"),verbose=TRUE){
#   require(parallel)
    if(verbose){
     cat(paste("Engaging warp drive using",type,"...\n",sep=" "))
    }
#
#   Start Cluster
#
    ### Set up the cluster
    cl <- parallel::makeCluster(parallel,type=type)
    ### initialize parallel random number streams
    if(is.null(seed)){
     parallel::clusterSetRNGStream(cl)
    }else{
     parallel::clusterSetRNGStream(cl,iseed=seed)
    }
    ### start each virtual machine with libraries loaded
    for(pkg in packagenames){
      attached <- parallel::clusterCall(cl, require, package=pkg, character.only=TRUE)      
    }
#
#   Run the jobs with Rmpi
#
    ### make sure that R has printed out console messages before go parallel
    utils::flush.console()
    return(cl)
}
#' @keywords internal
endparallel<-function(cl, type="MPI", finalize=TRUE, verbose=TRUE){
    parallel::stopCluster(cl)
# Activate the next line for Rmpi
#   if(finalize & type=="MPI"){Rmpi::mpi.finalize()}
    invisible()
}
