beginparallel<-function(parallel=1, type=NULL, seed=NULL, packagenames=c("size"),verbose=TRUE){
    ### parallel is wrapper for MPI or PVM (mosix only has PVM)
#   require(snow)
    require(parallel)
#
#   Start PVM if necessary
#
#   setDefaultClusterOptions(type="PVM")
#   setDefaultClusterOptions(type="MPI")
#   if(snow::getClusterOption("type")=="PVM"){
    if(is.null(type)){
     silentwarnings <- capture.output(try.rpvm<-require(rpvm, quietly=TRUE, warn.conflicts = FALSE))
     if(try.rpvm){
      snow::setDefaultClusterOptions(type="PVM")
      type <- "PVM"
      if(verbose){
       cat("Default warp drive is PVM ...\n")
      }
     }else{
      snow::setDefaultClusterOptions(type="MPI")
      type <- "MPI"
      if(verbose){
       cat("Default warp drive is MPI ...\n")
      }
     }
    }
    if(type=="PVM"){
     silentwarnings <- capture.output(try.rpvm<-require(rpvm, quietly=TRUE, warn.conflicts = FALSE))
     if(try.rpvm){
      if(verbose){
       cat("Engaging warp drive using PVM ...\n")
      }
#     require(rpvm)
      PVM.running <- try(rpvm::.PVM.config(), silent = TRUE)
      if(inherits(PVM.running, "try-error")){
       hostfile <- paste(Sys.getenv("HOME"), "/.xpvm_hosts", sep = "")
       if(file.exists(hostfile)){
        rpvm::.PVM.start.pvmd(hostfile)
       }else{
        rpvm::.PVM.start.pvmd()
       }
       cat("no problem... PVM started...\n")
      }
     }else{
      type <- "MPI"
      if(verbose){
       cat("PVM is not available. Engaging warp drive using MPI ...\n")
      }
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
    ### start each virtual machine with libraries loaded
    for(pkg in packagenames){
      attached <- clusterCall(cl, require, package=pkg, character.only=TRUE)      
    }
#
#   Run the jobs with rpvm or Rmpi
#
    ### make sure that R has printed out console messages before go parallel
    flush.console()
    return(cl)
}
endparallel<-function(cl, type=NULL, finalize=FALSE, verbose=TRUE){
    ### stop cluster and PVM (in case PVM is flakey)
    stopCluster(cl)
    if(snow::getClusterOption("type")=="PVM"){rpvm::.PVM.exit()}
    if(finalize & snow::getClusterOption("type")=="MPI"){Rmpi::mpi.finalize()}
#   if(type=="PVM"){
#    if(require("rpvm",character.only = TRUE)){
#     rpvm::.PVM.exit()
#    }else{
#     type <- "MPI"
#    }
#   }
#   if(type=="MPI"){
#    if(require("Rmpi",character.only = TRUE)){
#     Rmpi::mpi.finalize()
#    }
#   }
    invisible()
}
