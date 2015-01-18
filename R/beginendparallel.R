beginparallel<-function(parallel=1, type="MPI", seed=NULL, packagenames=c("size"),verbose=TRUE){
    require(parallel)
    if(verbose){
     cat("Engaging warp drive using MPI ...\n")
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
    flush.console()
    return(cl)
}
endparallel<-function(cl, type="MPI", finalize=TRUE, verbose=TRUE){
    parallel::stopCluster(cl)
    if(finalize & type=="MPI"){Rmpi::mpi.finalize()}
    invisible()
}
