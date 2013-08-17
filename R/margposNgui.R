# Author: JaneAc
#
# TODO/NOTES:
# margposN function is in the size.mle.R file
# should method and maxit be set by the user?
# I ordered the data based on waves
# NAs cause problems in the function. should they be removed as part of the function?
# Make output more useful: 1) return pmf as a variable than can be plotted? - Right now it's the qprob vector if return.all=TRUE)
	#add option to plot pmf?
# n=tabulate(s,nbin=K), error as with margpossize
###############################################################################

.makeMarginalPosteriorDistribution <- function() {
	
	JPanel <- J("javax.swing.JPanel")
	JButton <- J("javax.swing.JButton")
	JSeparator <- J("javax.swing.JSeparator")
	Dimension <- J("java.awt.Dimension")
	AnchorLayout <- J("org.rosuda.JGR.layout.AnchorLayout")
	SingletonDJList <- J("org.rosuda.deducer.toolkit.SingletonDJList")	
	ActionListener <- J("org.rosuda.deducer.widgets.event.RActionListener")		
	
	dialog <- new(Deducer::SimpleRDialog)
	dialog$setSize(400L, 400L)
	dialog$setTitle("Calculate Marginal Posterior Distribution of Population Size")
	
	varSel <- new(Deducer::VariableSelectorWidget)
	varSel$setRDataFilter("is.rds.data.frame")
	
	degreevar <- new(Deducer::SingleVariableWidget,"Degree Var",varSel)
	wavevar <- new(Deducer::SingleVariableWidget,"Wave Var",varSel) #wave or order?
	
	N_count <- new(Deducer::TextFieldWidget, "Pop. Size")
	N_count$setLowerBound(1)
	
	#K:
	max_s <- new(Deducer::TextFieldWidget, "Max Degree")
	max_s$setLowerBound(1)
	
	samples <- new(Deducer::TextFieldWidget, "# Samples")
	samples$setDefaultModel("100000")
	samples$setInteger(TRUE)
	samples$setLowerBound(1)

	#top
	addComponent(dialog, varSel, 50, 450, 550, 50)
	addComponent(dialog, degreevar, 50, 950, 200, 500)
	addComponent(dialog, wavevar, 250, 950, 400, 500, topType="REL",bottomType="REL")	
	
	#bottom
	addComponent(dialog, samples, 650, 325, 850, 75, topType="REL",bottomType="NONE")
	addComponent(dialog, N_count, 650, 625, 850, 375, topType="REL",bottomType="NONE")
	addComponent(dialog, max_s, 650, 925, 850, 675, topType="REL",bottomType="NONE")
	
	
	checkFunc <- function(x) {
		
		if (degreevar$getRModel() == "c()") 
			return("Please enter degree variable")
		if (wavevar$getRModel() == "") 
			return("Please enter wave variable")
		if (N_count$getModel()=="")
			return("Please enter population size")
		if (samples=="")
			return("Please enter number of samples")
		else("")
	
	}
	
	dialog$setCheckFunction(toJava(checkFunc))	
	
	
	runFunc <- function(x){
		
		"%+%" <- function(x, y) paste(x, y, sep = "")
		
		N <- N_count$getModel()
		
		deg<- unlist(strsplit(degreevar$getRModel(), "[\"]"))[2]
		
		wave<- unlist(strsplit(wavevar$getRModel(), "[\"]"))[2]		
		
		sbase1 <- varSel$getModel() %+% "$" %+% deg 
		sbase2 <- varSel$getModel() %+% "$" %+% wave
		s <- sbase1 %+% "[order(" %+% sbase2 %+% ")]" #leaving NAs for now [!is.na(" %+% sbase1 %+% "[order(" %+% sbase2 %+% ")])]"
		
		K <- max_s$getModel()
		if (K=="") {K = "max(s)"} 
		
		M <- samples$getModel() 
		
		method = "\"Nelder-Mead\"" #may update later to be set in dialog
		maxit = "500"
		if (method != "\"Nelder-Mead\"") {maxit = "100"}
		
		cmd <- "marg.pos.n<-margposN(" %+%
			" N = " %+% N %+%
			", s = " %+% s %+%
			", K = " %+% K %+%
			", M = " %+% M %+%
			", maxit = " %+% maxit %+%
			", method = " %+% method %+% 
			", parallel=1, seed=NULL,verbose=FALSE" %+%
			# n=tabulate(s,nbin=K), 
			", temp = 10, trace = 0, return.all=FALSE"  %+%
			")"
	
		#not entered through dialog:
		#parallel=1, seed, verbose, n, maxit, temp, trace, return.all
	    # prob = rep(1/K, K) , qprob = NULL are set internally
	
	execute(cmd)

	}	
	dialog$setRunFunction(toJava(runFunc))
	dialog
	}	