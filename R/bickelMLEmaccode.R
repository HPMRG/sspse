llhoodinC<-function(pop,s){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    N<-tabulate(pop)
    if(min(s)<1)print("Error - s values must be positive")
    n<-tabulate(s) #n is the vector of observed frequencies of each 
    if(length(n)>length(N)){print("Error - observed unit outside range")}
    Cret <- .C("bnw_llik",
              K=as.integer(length(N)),
              n=as.integer(length(s)),
              s=as.integer(s),
              snk=as.integer(n),
              Nk=as.double(N),
              llik=as.double(0), PACKAGE="size")
    Cret$llik
}

unposf<-function(mu,rho,N,s,n=tabulate(s)){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    if(length(n)>length(N)){print("Error - observed unit outside range")}
    Cret <- .C("bnw_unpos",
              K=as.integer(length(N)),
              n=as.integer(length(s)),
              s=as.integer(s),
              snk=as.integer(n),
              Nk=as.double(N),
              mu=as.double(mu),
              rho=as.double(rho),
              unpos=as.double(0), PACKAGE="size")
    Cret$unpos
}

dwarC<-function(N,mu,rho){
    # Waring PMF in C
    Cret <- .C("dwarC",
              N=as.integer(N),
              mu=as.double(mu),
              rho=as.double(rho),
              pmf=as.double(0), PACKAGE="size")
    Cret$pmf
}

llhoodfR<-function(N,s,n=tabulate(s)){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    n<-c(n,rep(0,times=length(N)-length(n)))
    sizes<-1:length(N)
    if(length(n)>length(N)){print("Error - observed unit outside range")}
    total_size<-sum(N*sizes)    
    obs_size<-cumsum(s) #cumulative sizes observed
    obs_size2<-c(0,obs_size)[1:length(obs_size)]
    result<-sum(lfactorial(N)-lfactorial(N-n))+sum(log(s/(total_size-obs_size2)))
    result
}

llhoodinCs<-function(pop,s){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    N<-tabulate(pop)
    if(min(s)<1)print("Error - s values must be positive")
#   if(length(n)>length(N)){print("Error - observed unit outside range")}
    Cret <- .C("bnw_s_llik",
              K=as.integer(length(N)),
              n=as.integer(length(s)),
              s=as.integer(s),
              Nk=as.double(N),
              llik=as.double(0), PACKAGE="size")
    Cret$llik
}



llhood<-function(pop,s){
    #this function takes a vector of population sizes and a vector s of 
    #sequential sizes of sampled units and returns a log likelihood value
    #s values must all be positive integers
    N<-tabulate(pop)
    if(min(s)<1)print("Error - s values must be positive")
    n<-tabulate(s) #n is the vector of observed frequencies of each 
    n<-c(n,rep(0,times=length(N)-length(n)))
    sizes<-1:length(N)
    if(length(n)>length(N)){print("Error - observed unit outside range")}
    total_size<-sum(N*sizes)    
    obs_size<-cumsum(s) #cumulative sizes observed
    obs_size2<-c(0,obs_size)[1:length(obs_size)]
    result<-sum(lfactorial(N)-lfactorial(N-n))+sum(log(s/(total_size-obs_size2)))
    result
}

llhoodfeed<-function(counts,sizes,s){
  # a function to feed values to llhood - more convenient input form
  #ord<-rev(order(sizes))
  #cts<-counts[ord]
  #ssiz<-sizes[ord]
  pop<-rep(sizes,times=counts)   
  llhood(pop,s)
  }

#brute force solving function.  requires min and max and does grid search - 2- dimensions only
bfsolve<-function(min,max,fn,...){#note: maximizes, also, exactly 2 dimensions
  dims<-length(min)
  best<-min
  bestscore<-fn(best,...)
  thisone<-min
  thisdim<-1
  for(j in min[thisdim]:max[thisdim]){
    thisdim<-1
    thisone[thisdim]<-j
    thisdim<-2
    for(i in min[thisdim]:max[thisdim]){
      thisone[thisdim]<-i
      #print(thisone)
      temp<-fn(thisone,...)
      if(temp>bestscore){best<-thisone; bestscore<-temp}
    }
}
print(bestscore)
best
}

# this function is a bit better - still brute force, but now loops over 
# arbitrary number of dimensions
bbfsolve<-function(minn,maxx,fn,...){#note: maximizes
  dims<-length(minn)
  best<-minn
  bestscore<-fn(best,...)
  thisone<-minn
  level<-1
  #print(thisone)
  wombat<-c(bestscore,best)
  #print(wombat)
  wombat<-looplower(level,minn,maxx,dims,wombat,fn,thisone,...) 
#  looplower(level,minn,maxx,dims,fn,best,thisone,bestscore,...) 
  bestscore<-wombat[1]
  best<-wombat[-1]
print(bestscore)
best
}

#subfunction.  note australia reference
looplower<-function(level,minn,maxx,dims,wombat,fn,thisone,...){
	#print(wombat)
	bestscore<-wombat[1]
	best<-wombat[-1]
#looplower<-function(level,minn,maxx,dims,fn,best,thisone,bestscore,...){
	#print(level)
	#print(minn)
	#print(maxx)
	#print(stuff)
	#print(dims)
	#print(thisone)
	for(i in minn[level]:maxx[level]){
		#print(thisone)
		thisone[level]<-i
		#print(thisone)
#		if(level<dims){looplower(level+1,minn,maxx,dims,fn,best,thisone,bestscore,...)
		if(level<dims){wombat<-looplower(level+1,minn,maxx,dims,wombat,fn,thisone,...)
			}else{temp<-fn(thisone,...); #print(temp)
      if(temp>bestscore){best<-thisone; bestscore<-temp; wombat<-c(bestscore,best); print(wombat)}  }
		}
		level<-level-1
	
	#print(wombat)
	wombat	
	}
