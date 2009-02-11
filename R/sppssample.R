llspps<-function(counts, sizes, s){
  #note counts and sizes refer to the counts of each size in the full
  # population.
  # s is the sample.
  if(min(s)<1)print("Error - s values must be positive")
  totalsize<-sum(counts*sizes)
  obs_size<-cumsum(s) #cumulative sizes observed
  obs_size2<-c(0,obs_size)[1:length(obs_size)]
  obscounts<-tabulate(s)[sizes]
  result<-sum(lfactorial(counts)-lfactorial(counts-obscounts))+sum(log(s/(totalsize-obs_size2)))
  result
}


lllspps<-function(lcountdif,sizes,s){
  #computes the log likelihood.  This is used to feed to optim.  The first
  # argument is the log of the
  # difference between the population counts and the sample counts.  This
  # formulation forces this difference
  # to be non-negative.
  countdif<-exp(lcountdif)
  counts<-countdif+tabulate(s)[sizes]
  llspps(counts,sizes,s)
}

sppssample<-function(counts,sizes,n){
   rstuff<-counts
   ss<-rep(0,n)
   for(i in 1:n){
     ss[i] <- sample(sizes, size=1,prob=rstuff*(sizes))
     rstuff[which(sizes==ss[i])] <- rstuff[which(sizes==ss[i])] - 1
   }
   ss
}

roundstoc<-function(vec){# takes a vector and makes it integers, keeping the total the same
        #vec<-c(1.2,1.3,2.5)
        target<-sum(vec)
        temp<-floor(vec)
        needed<-target-sum(temp)
        away<-vec-temp
        while(needed>.5){
                toget<-sample(c(1:length(vec)),size=1,prob=away)
                temp[toget]<-temp[toget]+1
                away<-vec-temp
                away[away<0]<-0
                needed<-needed-1
                }
        temp
        }
