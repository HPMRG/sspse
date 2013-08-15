# TODO/NOTES:
# should method and maxit be set by the user?
# I ordered the data based on waves and removed NAs
# ERROR: printing the following line repeatedly after calculations
# N=10000 minvalid=5000 nvalidhtn=0 nvalidhtd=1000000
# Author: JaneAc
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
		
		#not yet incorporated in function
		wave<- unlist(strsplit(wavevar$getRModel(), "[\"]"))[2]		
		
		sbase1 <- varSel$getModel() %+% "$" %+% deg 
		sbase2 <- varSel$getModel() %+% "$" %+% wave
		s <- sbase1 %+% "[order(" %+% sbase2 %+% ")][!is.na(" %+% sbase1 %+% "[order(" %+% sbase2 %+% ")])]"
		
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
			", prob = rep(1/K, K), parallel=1, seed=NULL,verbose=FALSE" %+%
			# n=tabulate(s,nbin=K), 
			", qprob = NULL, temp = 10, trace = 0, return.all=FALSE"  %+%
			")"
	
		#not entered through dialog: 
		#parallel=1, seed, verbose, n, qprob, maxit, temp, trace, return.all
	
	execute(cmd)

	}	
	dialog$setRunFunction(toJava(runFunc))
	dialog
	}	