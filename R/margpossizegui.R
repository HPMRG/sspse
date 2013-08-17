# Author: JaneAc
#
# NOTES:
# should N be set by user? currently it's not.
# a lot of output is just printed. should it be saved to a var instead (what's in result.all = True)?
# error: n=tabulate(s,nbin=K) as with margposN and posteriorsize
# Currently there's error: "Error: The population counts should not be exceeded by the sample counts."
###############################################################################


.makeMarginalPosteriorSize <- function() {
	
	dialog <- new(Deducer::SimpleRDialog)
	dialog$setSize(400L, 400L)
	dialog$setTitle("Calculate Marginal Posterior Size Distribution")
	
	varSel <- new(Deducer::VariableSelectorWidget)
	varSel$setRDataFilter("is.rds.data.frame")
	
	degreevar <- new(Deducer::SingleVariableWidget,"Degree Var",varSel)
	wavevar <- new(Deducer::SingleVariableWidget,"Wave Var",varSel) #wave or order?
	
	#K:
	max_s <- new(Deducer::TextFieldWidget, "Max Degree")
	max_s$setLowerBound(1)
	
	#M:
	samples <- new(Deducer::TextFieldWidget, "# Samples")
	samples$setDefaultModel("100000")
	samples$setInteger(TRUE)
	samples$setLowerBound(1)
	
	#top
	addComponent(dialog, varSel, 50, 450, 550, 50)
	addComponent(dialog, degreevar, 50, 950, 200, 500)
	addComponent(dialog, wavevar, 250, 950, 400, 500, topType="REL",bottomType="REL")	
	
	#bottom
	addComponent(dialog, samples, 650, 450, 850, 75, topType="REL",bottomType="NONE")
	addComponent(dialog, max_s, 650, 925, 850, 550, topType="REL",bottomType="NONE")
	
	
	checkFunc <- function(x) {
		
		if (degreevar$getRModel() == "c()") 
			return("Please enter degree variable")
		if (wavevar$getRModel() == "c()") 
			return("Please enter wave variable")
		if (samples=="")
			return("Please enter number of samples")
		else("")
		
	}
	dialog$setCheckFunction(toJava(checkFunc))	
	
	runFunc <- function(x){
		
		"%+%" <- function(x, y) paste(x, y, sep = "")
	
		deg<- unlist(strsplit(degreevar$getRModel(), "[\"]"))[2]
		
		wave<- unlist(strsplit(wavevar$getRModel(), "[\"]"))[2]		
		
		sbase1 <- varSel$getModel() %+% "$" %+% deg 
		sbase2 <- varSel$getModel() %+% "$" %+% wave
		s <- sbase1 %+% "[order(" %+% sbase2 %+% ")]" #leaving NAs for now [!is.na(" %+% sbase1 %+% "[order(" %+% sbase2 %+% ")])]"
		
		
		K <- max_s$getModel()
		if (K=="") {K = "max(s)"}
		
		M <- samples$getModel() 
		
		N = "trunc(length(s)*seq(1.1,4,length=10)+1)" #not currently set by user
		
		cmd <- "marg.pos.size<-margpossize(" %+%	
		"s = " %+% s %+%
		", N = " %+% N %+%
		", K = " %+% K %+%
		", M = " %+% M %+%
		", parallel=1,seed=NULL, verbose=FALSE" %+% #error: n=tabulate(s,nbin=K), 
		", return.all=FALSE)"
		
		
	execute(cmd)
			

	}	
	dialog$setRunFunction(toJava(runFunc))
	dialog
	
	
}
