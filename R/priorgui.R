# TODO: Add comment
# 
# Author: JaneAc
###############################################################################


.makePriorDistribution <- function() {

	#column1	
	dialog <- new(Deducer::SimpleRDialog)
	dialog$setSize(500L, 500L)
	dialog$setTitle("Calculate Prior Distribution of Population Size")
	
	samplesize <- new(Deducer::TextFieldWidget, "Sample Size")
	samplesize$setInteger(TRUE)
	samplesize$setLowerBound(1)
	
	maxN <- new(Deducer::TextFieldWidget, "Population Max")
	maxN$setInteger(TRUE)
	maxN$setLowerBound(1)
	
	typedist <-new(Deducer::ButtonGroupWidget, "Distribution Type", c("proportion","nbinom","continuous","pln","flat"))
	typedist$setToolTipText("nbinom = Negative Binomial; pln = Poison-log-normal")
	
	plotbox <- new(Deducer::CheckBoxesWidget,.jarray("Plot Distribution"))
	setSize(plotbox,400,75)	
	#plotbox$setDefaultModel(c("Plot Distribution"))
	
	#column2
	priormean <- new(Deducer::TextFieldWidget, "Prior Mean")
	priormean$setLowerBound(1)
	
	priorsd <- new(Deducer::TextFieldWidget, "Prior Standard Deviation")
	priorsd$setLowerBound(0)
	
	priormed <- new(Deducer::TextFieldWidget, "Prior Median")
	priormed$setLowerBound(1)
	
	priormode <- new(Deducer::TextFieldWidget, "Prior Mode")
	priormed$setLowerBound(1)

	priorquartiles <- new(Deducer::TextFieldWidget, "Prior Quartiles (25%, 75%)")
	#set to require acoordinate pair
	
	#ignore for now - user can set on command line
	#prioralpha <- new(Deducer::TextFieldWidget, "Alpha")
	#effectivedf <- new(Deducer::TextFieldWidget, "Effective Prior DF")
	#effectivedf$setDefaultModel("1")
	
	priormodeprop <- new(Deducer::TextFieldWidget, "Proportion Prior Mode")
	priormodeprop$setLowerBound(0)
	priormodeprop$setUpperBound(1)
	priormodeprop$setDefaultModel(".5")
	
	priormedprop <- new(Deducer::TextFieldWidget, "Proportion Prior Median")
	priormedprop$setLowerBound(0)
	priormedprop$setUpperBound(1)
	
	#left column 1
	addComponent(dialog, samplesize, 50, 350, 150, 100, bottomType = "NONE")
    addComponent(dialog, maxN, 175, 350, 275, 100, bottomType = "NONE")
	addComponent(dialog, typedist, 325, 450, 675, 100, bottomType = "NONE")
	addComponent(dialog, plotbox, 675, 450, 750, 100, bottomType = "NONE")

	#right column 2
	addComponent(dialog, priormean, 50, 900, 150, 500, bottomType = "NONE")
	addComponent(dialog, priorsd, 175, 900, 275, 500, bottomType = "NONE")
	addComponent(dialog, priormed, 300, 900, 400, 500, bottomType = "NONE")
	addComponent(dialog, priormode, 425, 900, 525, 500, bottomType = "NONE")
	addComponent(dialog, priorquartiles, 550, 900, 650, 500, bottomType = "NONE")
	addComponent(dialog, priormodeprop, 675, 900, 775, 500, bottomType = "NONE")
	addComponent(dialog, priormedprop, 800, 900, 900, 500, bottomType = "NONE")
	
	checkFunc <- function(x) {
		if (samplesize$getModel() == "") 
			return("Please enter sample size")
		if (maxN$getModel() == "") 
			return("Please enter population max")
		if (priorquartiles$getModel()!="" && !length(strsplit(priorquartiles$getModel(),",")[[1]])%in%c(0,2))
			#return("Quartile entry should be empty or vector of length two")
			return(strsplit(priorquartiles$getModel(),",")[[1]])
		if (priormodeprop$getModel()>1)
			return("Prior propotion mode must be between 0 and 1")
		if (priormedprop$getModel()>1)
			return("Prior propotion median must be between 0 and 1")
		else("")

		#if type is flat or pln, grey out everything but n and maxN
		#if type is proportioon, grey out mode.prior.size
		
		#tmp <- weightPanel$check()
		#if (tmp != "") 
		#	return(tmp)
		#if (sub$getModel() != "") {
		#	data <- varSel$getModel()
		#	subExp <- sub$getModel()
		#	if (!J("RDSAnalyst.SubsetWidget")$isValidSubsetExp(subExp, 
		#			data)) 
		#		return("Invalid subset")
		
}

	dialog$setCheckFunction(toJava(checkFunc))	

	runFunc <- function(x){
	
		"%+%" <- function(x, y) paste(x, y, sep = "")
		n <- samplesize$getModel()
		max_N <- maxN$getModel()
		dist_type <- switch(typedist$getModel(),
				"proportion"="proportion",  #Decided to leave names as is, but can update this later
				"nbinom"="nbinom",
				"continuous"="continuous",
				"pln"="pln",
				"flat" = "flat" )
		prior.mean <- "NULL"
			if (priormean$getModel()!="") {prior.mean = priormean$getModel()}
		prior.SD <- "NULL"
			if (priorsd$getModel()!="") {prior.SD= priorsd$getModel()}
		prior.med <- "NULL"
			if (priormed$getModel()!="") {prior.med = priormed$getModel()}
		prior.mode <- "NULL"
			if (priormode$getModel()!="") {prior.mode = priormode$getModel()}
		prior.quart <- "NULL"
		quarts <- "c("%+%priorquartiles$getModel()%+%")"
		#add hovertext to ensure proper entry form and/or add check for extra parenths
			if (priorquartiles$getModel()!="") {prior.quart = quarts}
		prior.prop.mode <- "NULL"
			if (priormodeprop$getModel()!="") {prior.prop.mode = priormodeprop$getModel()}
		prior.prop.med <- "NULL"
			if (priormedprop$getModel()!="") {prior.prop.med = priormedprop$getModel()}	
			
		cmd <- "dsp <- dsizeprior(" %+% n %+% ", type=\"" %+% dist_type %+% "\", mean.prior.size =" %+% prior.mean %+%
		", maxN=" %+% max_N %+%	", sd.prior.size=" %+% 
		prior.SD %+% ", mode.prior.sample.proportion=" %+% prior.prop.mode %+% 
		", median.prior.sample.proportion=" %+% prior.prop.med %+% 
		", median.prior.size=" %+% prior.med %+% ", mode.prior.size =" %+% 
		prior.mode %+% ", quartiles.prior.size=" %+% prior.quart %+% 
		", effective.prior.df=1, alpha = NULL, beta = NULL, log = FALSE, maxbeta = 100, maxNmax = 2e+05, verbose = TRUE" 
	
		
		
		if(plotbox$getModel()$size()>0) { #return info and plot instead of pmf vectors
			cmd <- cmd %+% "); dsp[3:14]; plot(dsp$x,dsp$lprior,main = \"Prior Pop. Size Distribution (" %+% dist_type %+% ")\")"
			}
			else (cmd <- cmd %+% ");dsp\n")
		execute(cmd)
	}	
	dialog$setRunFunction(toJava(runFunc))
	dialog
}