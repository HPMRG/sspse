# Author: JaneAc

#TO DO:
#na.rm's in quantile func?
#optim(par bug
#samping defaults should be different for posterior size
#mention (in tooltip?) that prior data can be calculated via prior distribution dialog
#add info buttons on subdialogs with screenshots of annotated versions of the dialogs


###############################################################################

###############################################################################

.makePosteriorDistribution <- function() {
	
	JPanel <- J("javax.swing.JPanel")
	JButton <- J("javax.swing.JButton")
	JSeparator <- J("javax.swing.JSeparator")
	Dimension <- J("java.awt.Dimension")
	AnchorLayout <- J("org.rosuda.JGR.layout.AnchorLayout")
	SingletonDJList <- J("org.rosuda.deducer.toolkit.SingletonDJList")	
	ActionListener <- J("org.rosuda.deducer.widgets.event.RActionListener")	
	MouseListener <- J("org.rosuda.deducer.widgets.event.RMouseListener")	
	
	
	#Components
	#Buttons
	#Arrange Components
	
	
	#Components
	
	##TOP: choose degree and optional status variables. 
	dialog <- new(Deducer::SimpleRDialog)
	dialog$setSize(500L, 550L)
	dialog$setTitle("Calculate Posterior Distribution of Population Size")
	
	varSel <- new(Deducer::VariableSelectorWidget)
	varSel$setRDataFilter("is.rds.data.frame")
	
	degreevar <- new(Deducer::SingleVariableWidget,"Degree Var",varSel)
	diseasevar <- new(Deducer::SingleVariableWidget,"Status Var",varSel)
	wavevar <- new(Deducer::SingleVariableWidget,"Wave Var",varSel) #wave or order?
	
	diseasebox <- new(Deducer::CheckBoxesWidget,.jarray("Use Status Variable")) 	##choose to use posteriorsize or posteriordisease function via checkbox (diseasebox)
	
	
	##SAMPLING
	intervalsize <- new(Deducer::TextFieldWidget, "Interval")
	intervalsize$setDefaultModel("10")
	intervalsize$setLowerBound(1)
	
	burn <- new(Deducer::TextFieldWidget, "Burnin")
	burn$setDefaultModel("5000")
	burn$setLowerBound(0)
	
	samples <- new(Deducer::TextFieldWidget, "Number of Samples")
	samples$setDefaultModel("1000")
	samples$setLowerBound(1)
	
	
	#POP PRIOR:	
	max_N <- new(Deducer::TextFieldWidget, "Pop. Max")
	max_N$setLowerBound(1)

	priormed <- new(Deducer::TextFieldWidget, "Median")
	priormed$setLowerBound(1)
	
	priormean <- new(Deducer::TextFieldWidget, "Mean")
	priormean$setLowerBound(1)
	
	priorsd <- new(Deducer::TextFieldWidget, "SD")
	priorsd$setLowerBound(0)
	
	priormode <- new(Deducer::TextFieldWidget, "Mode")
	priormode$setLowerBound(1)
	
	quarts <- new(Deducer::TextFieldWidget, "Quartiles")
	
	priormodeprop <- new(Deducer::TextFieldWidget, "Prop. Mode")
	priormodeprop$setLowerBound(0)
	priormodeprop$setUpperBound(1)
	priormodeprop$setDefaultModel(".5")
	
	types = c("Proportion","Flat","Neg-binom","Poisson-log-norm")
	typedist <-new(Deducer::ComboBoxWidget, types)
	typedist$setTitle("Prior Dist. Type", TRUE)
	typedist$setDefaultModel("proportion")
	
	
	##DEGREE
	maxDeg <- new(Deducer::TextFieldWidget, "Max Degree (K)") #= K, default = round(quantile(s,0.80))
	maxDeg$setLowerBound(1)
	#maxDegListen <- new(MouseListener)
	#maxDeg$setToolTipText("The maximum degree for an individual. This is usually calculated as twice the maximum observed degree.")
	#maxDeg$addMouseListener(maxDegListen)
	
	degs1 = c("Conway-Maxwell-Poisson","Neg-binom","Poisson-log-normal")
	degreedist <-new(Deducer::ComboBoxWidget, degs1)
	degreedist$setTitle("Degree Distribution Type", TRUE)
	degreedist$setDefaultModel("Conway-Maxwell-Poisson")
	
	priordegreemean <- new(Deducer::TextFieldWidget, "Mean Degree")
	priordegreemean$setLowerBound(0)
	
	priordegreeSD <- new(Deducer::TextFieldWidget, "Std. Dev.")
	priordegreeSD$setLowerBound(0)
	
	priordegreemean0 <- new(Deducer::TextFieldWidget, "Mean, Status = 0")
	priordegreemean0$setLowerBound(0)
	priordegreemean0$setDefaultModel("7")
	
	priordegreemean1 <- new(Deducer::TextFieldWidget, "Mean, Status = 1")
	priordegreemean1$setLowerBound(0)
	priordegreemean1$setDefaultModel("7")
	
	dispers <- new(Deducer::TextFieldWidget, "Dispersion")
	dispers$setLowerBound(0)
	
	
	#Buttons 
	mcbutton <- new(JButton,"Sampling Prefs")
	mcbutton$setToolTipText("Set burnin period, sampling interval and number of samples for MCMC sampling.")
	mcsubdialog <- new(SimpleRSubDialog,dialog,"Set MCMC Sampling Parameters")
	setSize(mcsubdialog,300L,200L)
	
	mcactionFunction <- function(cmd,ActionEvent){
		mcsubdialog$setLocationRelativeTo(mcbutton)
		mcsubdialog$run()
	}
	mclistener <- new(ActionListener)
	mclistener$setFunction(toJava(mcactionFunction))
	mcbutton$addActionListener(mclistener)
	
	degbutton <- new(JButton,"Prior Degree Data")
	degbutton$setToolTipText("Only use Status = 0 and Status = 1 fields if \"Use Status Var\" is checked.")	
	degsubdialog <- new(SimpleRSubDialog,dialog,"Prior Degree Data")
	setSize(degsubdialog,300L,400L)
	
	degactionFunction <- function(cmd,ActionEvent){
		degsubdialog$setLocationRelativeTo(degbutton)
		degsubdialog$run()
	}
	deglistener <- new(ActionListener)
	deglistener$setFunction(toJava(degactionFunction))
	degbutton$addActionListener(deglistener)
	
	popbutton <- new(JButton,"Prior Population Data")
	
	popsubdialog <- new(SimpleRSubDialog,dialog,"Prior Population Data")
	setSize(popsubdialog,300L,400L)
	
	popactionFunction <- function(cmd,ActionEvent){
		popsubdialog$setLocationRelativeTo(popbutton)
		popsubdialog$run()
	}
	poplistener <- new(ActionListener)
	poplistener$setFunction(toJava(popactionFunction))
	popbutton$addActionListener(poplistener)
	
	
	#Arrange Components:

	#top 
	addComponent(dialog, varSel, 50, 450, 550, 50)
	#variable boxes
	addComponent(dialog, degreevar, 50, 950, 175, 500, topType="ABS",bottomType="NONE")
	#setSize(diseasevar,400,125)
	addComponent(dialog, diseasevar, 200, 950, 325, 500, topType="REL",bottomType="REL")
	#checkbox
	addComponent(dialog, diseasebox, 325, 950, 425, 625, bottomType = "REL")
	
	addComponent(dialog, wavevar, 425, 950, 550, 500, topType="REL",bottomType="REL")
	

	#buttons
	addComponent(dialog,mcbutton,625,450,725,50)
	addComponent(dialog,degbutton,775,450,875,50)
	addComponent(dialog,popbutton,625,950,725,550)	
	
	#Sampling button
		#replaced panel with subdialog (JC)
		#samplepanel <- new(JPanel)
		#samplepanel$setPreferredSize(new(Dimension,350L,250L))
		#samplepanel$setSize(new(Dimension,350L,250L))
		#samplepanel$setBorder(J("javax.swing.BorderFactory")$createTitledBorder("MCMC Sampling"))
		#samplepanel$setLayout(new(AnchorLayout))
		#addComponent(dialog, samplepanel, 375, 450, 625, 50,bottomType="REL")		
	
		#addComponent(samplepanel, burn, 100, 475, 450, 50, bottomType = "NONE")
		#addComponent(samplepanel, intervalsize, 100, 950, 450, 525, bottomType = "NONE")
		#addComponent(samplepanel, samples, 500, 950, 950, 50, bottomType = "NONE")
	

	
	addComponent(mcsubdialog, burn, 100, 475, 450, 50, bottomType = "NONE")
	addComponent(mcsubdialog, intervalsize, 100, 950, 450, 525, bottomType = "NONE")
	addComponent(mcsubdialog, samples, 500, 950, 850, 50, bottomType = "NONE")



	#degreeinfo
	degreepanel <- new(JPanel)
	degreepanel$setBorder(J("javax.swing.BorderFactory")$createTitledBorder("Prior Degree Data"))
	degreepanel$setLayout(new(AnchorLayout))
	#addComponent(dialog, degreepanel, 650, 450, 900, 50,bottomType="REL")		

	addComponent(degsubdialog, maxDeg, 100, 475, 200, 50, bottomType = "NONE")
	addComponent(degsubdialog, dispers, 100, 950, 200, 525, bottomType = "NONE")
	addComponent(degsubdialog, priordegreemean, 250, 475, 350, 50, bottomType = "NONE")
	addComponent(degsubdialog, priordegreeSD, 250, 950, 350, 525, bottomType = "NONE")
	addComponent(degsubdialog, priordegreemean0, 400, 475, 500, 50, bottomType = "NONE")
	addComponent(degsubdialog, priordegreemean1, 400, 950, 500, 525, bottomType = "NONE")
	addComponent(degsubdialog, degreedist, 575, 950, 800, 50, bottomType = "REL")
	
	
	#column 2
	
#prior data panel
		#replaced panel with subdialog
		#priorpanel <- new(JPanel)
		#priorpanel$setBorder(J("javax.swing.BorderFactory")$createTitledBorder("Prior Population Data"))
		#priorpanel$setLayout(new(AnchorLayout))
		#addComponent(dialog, priorpanel, 375, 950, 900, 500,bottomType="REL")		
	
	addComponent(popsubdialog, max_N, 50, 475, 200, 50, bottomType = "REL")
	addComponent(popsubdialog, priormed, 50, 950, 200, 525, bottomType = "REL")
	
	addComponent(popsubdialog, priormean, 250, 475, 400, 50, bottomType = "REL")
	addComponent(popsubdialog, priorsd, 250, 950, 400, 525, bottomType = "REL")
	
	addComponent(popsubdialog, priormode, 450, 475, 600, 50, bottomType = "REL")
	addComponent(popsubdialog, quarts, 450, 950, 600, 525, bottomType = "REL")
	
	addComponent(popsubdialog, priormodeprop, 650, 400, 800, 50, bottomType = "REL")
	addComponent(popsubdialog, typedist, 650, 950, 800, 450, bottomType = "REL")

	checkFunc <- function(x) {
		
		if (degreevar$getRModel() == "c()") 
			return("Please enter degree variable")
		if (wavevar$getRModel() == "c()") 
			return("Please enter wave variable")
		if (diseasevar$getRModel() == "c()" && diseasebox$getModel()$size()>0) 
			return("Please enter status variable")
		if (quarts$getModel()!="" && !length(strsplit(quarts$getModel(),",")[[1]])%in%c(0,2))
			return("Quartile entry should be empty or of form low,high e.g 2000,5000")
		if (priormodeprop$getModel()>1)
			return("Prior proportion mode must be between 0 and 1")
		#if (priormedprop$getModel()>1)
		#	return("Prior proportion median must be between 0 and 1")
		else return("")
		
	}
	dialog$setCheckFunction(toJava(checkFunc))	
	
	runFunc <- function(x){
			
		"%+%" <- function(x, y) paste(x, y, sep = "")
		
		deg<- unlist(strsplit(degreevar$getRModel(), "[\"]"))[2]
		
		#not yet incorporated in function
		wave<- unlist(strsplit(wavevar$getRModel(), "[\"]"))[2]		
		
		s <- varSel$getModel() %+% "$" %+% deg
		
		median.prior.size="NULL"
		if (priormed$getModel()!="") {median.prior.size = priormed$getModel()}
		interval <- "10"
		if (intervalsize$getModel()!="") {interval = intervalsize$getModel()}
		burnin <- "5000"
		if (burn$getModel()!="") {burnin = burn$getModel()}
		maxN = "NULL"
		
		if (max_N$getModel()!="") {maxN = max_N$getModel()}
		s2=eval(parse(text=s))
		K=as.character(round(quantile(s2,0.80,na.rm=TRUE))) #should this be double the max observed instead??
		if (maxDeg$getModel()!="") {K = maxDeg$getModel()}
		samplesize <- "1000"
		if (samples$getModel()!="") {samplesize = samples$getModel()}
		
		quartiles.prior.size <- "NULL"
		quarts2 <- quarts$getRModel()
		#add hovertext to ensure proper entry form and/or add check for extra parenths?
		if (quarts2!="\"\"") {quartiles.prior.size = quarts2}
		
		mean.prior.size <- "NULL"
		if (priormean$getModel()!="") {mean.prior.size = priormean$getModel()}
		mode.prior.size <- "NULL"
		if (priormode$getModel()!="") {mode.prior.size = priormode$getModel()}
		
		sd.prior.size <- "NULL"
		if (priorsd$getModel()!="") {sd.prior.size = priorsd$getModel()}
		mode.prior.sample.proportion <- ".5"
		if (priormodeprop$getModel()!="") {mode.prior.sample.proportion = priormodeprop$getModel()}
		mean.prior.degree <- "NULL"
		if (priordegreemean$getModel()!="") {mean.prior.degree = priordegreemean$getModel()}
		sd.prior.degree <- "NULL"
		if (priordegreeSD$getModel()!="") {sd.prior.degree = priordegreeSD$getModel()}
		dispersion <- "0"
		if (dispers$getModel()!="") {dispersion = dispers$getModel()}
		
		mean0.prior.degree = priordegreemean0$getModel()
		mean1.prior.degree = priordegreemean1$getModel()
		
		
		priorsizedistribution <- switch(typedist$getModel(),
				"Proportion"="proportion",
				"Neg-binom"="nbinom",
				"Poisson-log-norm"="pln",
				"Flat" = "flat" )
		
		priordegreedistribution <- switch(degreedist$getModel(),
				"Conway-Maxwell-Poisson" = "cmp", 
				"Neg-binom"="nbinom",
				"Poisson-log-norm"="pln")
		
		#no entry field
		alpha = "NULL"
		df.mean.prior="1"
		df.sd.prior = "5"
		Np = "0"
		#effective.prior.df=1
			#nk=tabulate(s,nbin=K),
		#n=length(s),
		#muproposal=0.1, 
		#sigmaproposal=0.15, 
		#burnintheta=500,
		#parallel=1, 
		#parallel.type="PVM",
		#seed=NULL
		#verbose=TRUE
		
		if (diseasebox$getModel()$size()<=0) {
		cmd <- "posize <- posteriorsize(s =" %+% s %+% #", median.prior.size=" %+% median.prior.size %+% ", interval =" %+% intervalsize %+%
					", burnin=" %+% burnin %+%
					", maxN=" %+% maxN %+% ", K=" %+% K %+% 
					", samplesize=" %+% samplesize %+% 
					", quartiles.prior.size=" %+% quartiles.prior.size %+% 
					", mean.prior.size =" %+% mean.prior.size %+% 
					", mode.prior.size=" %+% mode.prior.size %+% 
					", priorsizedistribution= \"" %+% priorsizedistribution %+% "\"" %+%
					", effective.prior.df= 1" %+% 
					", sd.prior.size =" %+%	sd.prior.size %+% 
					", mode.prior.sample.proportion=" %+% mode.prior.sample.proportion %+% 
					", alpha=" %+% alpha %+% 
					", degreedistribution= \"" %+% priordegreedistribution %+% "\"" %+%
					", mean.prior.degree=" %+% mean.prior.degree %+% 
					", sd.prior.degree=" %+% sd.prior.degree %+% 
					", df.mean.prior =" %+% df.mean.prior %+% 
					", df.sd.prior=" %+% df.sd.prior %+% 
					", Np=" %+% Np %+% 
					", dispersion=" %+% dispersion %+% 
				#", nk=tabulate(s,nbin=K)" %+%  This line causes crash when included excplicitly
					", n=length(s)" %+% 
					", muproposal=0.1, sigmaproposal=0.15, burnintheta=500" %+% 
					", parallel=1, parallel.type=\"PVM\", seed=NULL, verbose=TRUE" %+% 
					")"
			cmd <- cmd # %+% ";posize\n"
			
		}
		
		else { #use posteriordisease function				
				
				dis2<- unlist(strsplit(diseasevar$getRModel(), "[\"]"))[2]		
				dis <- varSel$getModel() %+% "$" %+% dis2
				Np0 = "0"
				Np1 = "0"
				burnintheta="500"	
				cmd <- "podisease <- posteriordisease(s = " %+% s %+% ", dis = " %+% dis %+%
						", mean0.prior.degree = " %+% mean0.prior.degree %+%
						", mean1.prior.degree = " %+% mean1.prior.degree %+%
						", sd.prior.degree = " %+% sd.prior.degree %+%
						", df.mean.prior = " %+% df.mean.prior %+% 
						", df.sd.prior=" %+% df.sd.prior %+% 
						", Np0= " %+% Np0 %+% 
						", Np1= " %+% Np1 %+% 
						", samplesize= " %+% samplesize %+%	
						", burnin= " %+% burnin %+%	
						", interval= " %+% interval %+%	
						", burnintheta= " %+% burnintheta %+%	
						", priorsizedistribution= \"" %+% priorsizedistribution %+% "\"" %+%
						", mean.prior.size = " %+% mean.prior.size %+% 
						", sd.prior.size = " %+%	sd.prior.size %+% 
						", mode.prior.sample.proportion = " %+% mode.prior.sample.proportion %+% 
						", median.prior.size = " %+% median.prior.size %+% 
						", mode.prior.size =" %+% mode.prior.size %+% 
						", quartiles.prior.size = " %+% quartiles.prior.size %+%
						", effective.prior.df = 1" %+% 
						", alpha = " %+% alpha %+% 
						", degreedistribution= \"" %+% priordegreedistribution %+% "\"" %+%
						", maxN = " %+% maxN %+% 
						", K = " %+% K %+% 
						", n = length(s)" %+%
						", dispersion = " %+% dispersion %+%
						", nk0 = tabulate(s[dis==0],nbin=K), nk1 = tabulate(s[dis==1],nbin = K)" %+%
						", muproposal = 0.1, sigmaproposal = 0.15, parallel = 1, seed = NULL" %+%
						", verbose = TRUE"  %+%
						")"

				cmd <- cmd # %+% ";podisease\n"
						}
		execute(cmd)

	}	
	dialog$setRunFunction(toJava(runFunc))
	dialog
}