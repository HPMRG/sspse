# TODO: Add comment
# 
# Author: JaneAc

#questions:
#na.rm's in quantile func?
#not working with nyjazz data, but ok with fauxmadronadata

###############################################################################

###############################################################################

.makePosteriorDistribution <- function() {
	
	JPanel <- J("javax.swing.JPanel")
	
	JSeparator <- J("javax.swing.JSeparator")
	Dimension <- J("java.awt.Dimension")
	AnchorLayout <- J("org.rosuda.JGR.layout.AnchorLayout")
	SingletonDJList <- J("org.rosuda.deducer.toolkit.SingletonDJList")	
	#RActionListener <- J("org.rosuda.deducer.widgets.event.RActionListener")
	
	
	#top and column 1 - choose degree and optional status variables. 
	#choose to use posteriorsize or posteriordisease function via checkbox (diseasebox)
	dialog <- new(Deducer::SimpleRDialog)
	dialog$setSize(500L, 700L)
	dialog$setTitle("Calculate Posterior Distribution of Population Size")
	
	varSel <- new(Deducer::VariableSelectorWidget)
	varSel$setRDataFilter("is.rds.data.frame")
	
	degreevar <- new(Deducer::SingleVariableWidget,"Degree Var",varSel)
	
	diseasevar <- new(Deducer::SingleVariableWidget,"Status Var",varSel)
	
	diseasebox <- new(Deducer::CheckBoxesWidget,.jarray("Use Status Variable"))
	
	#sampling information
	intervalsize <- new(Deducer::TextFieldWidget, "Interval")
	intervalsize$setDefaultModel("10")
	intervalsize$setLowerBound(1)
	
	burn <- new(Deducer::TextFieldWidget, "Burnin")
	burn$setDefaultModel("5000")
	burn$setLowerBound(0)
	
	samples <- new(Deducer::TextFieldWidget, "Number of Samples")
	samples$setDefaultModel("1000")
	samples$setLowerBound(1)
	
	
	#column 2
	#prior information:
	
	#can be calculated via prior distribution dialog
	
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
	
	priormodeprop <- new(Deducer::TextFieldWidget, "Proportion Mode")
	priormodeprop$setLowerBound(0)
	priormodeprop$setUpperBound(1)
	priormodeprop$setDefaultModel(".5")
	
	types = c("proportion","flat","nbinom","pln")
	typedist <-new(Deducer::ComboBoxWidget, types)
	typedist$setTitle("Prior Distribution Type", TRUE)
	#default automatically set to proportion
	
	#degree info
	maxDeg <- new(Deducer::TextFieldWidget, "Max Deg. (K)") #= K, default = round(quantile(s,0.80))
	maxDeg$setLowerBound(1)
	maxDeg$setToolTipText("The maximum degree for an individual. This is usually calculated as twice the maximum observed degree.")
	
	degreedist <-new(Deducer::ButtonGroupWidget, "Degree Distribution Type", c("cmp","nbinom","pln"))
	degreedist$setDefaultModel("CMP")
	#below isn't working. Doesn't work on ButtonGroupWidget?
	#degreedist$setToolTipText("cmp = Conway-Maxwell-Poisson; nbinom = Negative Binomial; pln = Poisson-log-normal")
	
	priordegreemean <- new(Deducer::TextFieldWidget, "Mean Degree")
	priordegreemean$setLowerBound(0)
	
	priordegreeSD <- new(Deducer::TextFieldWidget, "SD Degree")
	priordegreeSD$setLowerBound(0)
	
	dispers <- new(Deducer::TextFieldWidget, "Dispersion")
	dispers$setLowerBound(0)
	
	
	#plotbox <- new(Deducer::CheckBoxesWidget,.jarray("Plot Distribution"))
	#setSize(plotbox,400,75)	
	#plotbox$setDefaultModel(c("Plot Distribution"))
	
	#column2
	priormean <- new(Deducer::TextFieldWidget, "Mean")
	priormean$setLowerBound(1)

	
	priormode <- new(Deducer::TextFieldWidget, "Mode")
	priormed$setLowerBound(1)
	
	priorquartiles <- new(Deducer::TextFieldWidget, "Prior Quartiles (25%, 75%)")
	#set to require acoordinate pair
	
	#ignore for now - user can set on command line
	#prioralpha <- new(Deducer::TextFieldWidget, "Alpha")
	#effectivedf <- new(Deducer::TextFieldWidget, "Effective Prior DF")
	#effectivedf$setDefaultModel("1")
	
	#priormedprop <- new(Deducer::TextFieldWidget, "Prop. Med")
	#priormedprop$setLowerBound(0)
	#priormedprop$setUpperBound(1)
	
	#top and column 1
	addComponent(dialog, varSel, 50, 450, 350, 50)
	#variable boxes
	addComponent(dialog, degreevar, 50, 950, 125, 500, topType="ABS",bottomType="NONE")
	setSize(diseasevar,400,100)
	addComponent(dialog, diseasevar, 150, 950, 250, 500, topType="REL",bottomType="REL")
	#checkbox
	addComponent(dialog, diseasebox, 275, 950, 350, 600, bottomType = "REL")
	
	#column 1
	samplepanel <- new(JPanel)
	#samplepanel$setPreferredSize(new(Dimension,350L,250L))
	#samplepanel$setSize(new(Dimension,350L,250L))
	samplepanel$setBorder(J("javax.swing.BorderFactory")$createTitledBorder("MCMC Sampling"))
	samplepanel$setLayout(new(AnchorLayout))
	addComponent(dialog, samplepanel, 375, 450, 625, 50,bottomType="REL")		
	
	addComponent(samplepanel, burn, 100, 475, 450, 50, bottomType = "NONE")
	addComponent(samplepanel, intervalsize, 100, 950, 450, 525, bottomType = "NONE")
	addComponent(samplepanel, samples, 500, 950, 950, 50, bottomType = "NONE")
	
	#degreeinfo
	degreepanel <- new(JPanel)
	degreepanel$setBorder(J("javax.swing.BorderFactory")$createTitledBorder("Prior Degree Data"))
	degreepanel$setLayout(new(AnchorLayout))
	addComponent(dialog, degreepanel, 650, 450, 900, 50,bottomType="REL")		

	addComponent(degreepanel, maxDeg, 100, 475, 450, 50, bottomType = "NONE")
	addComponent(degreepanel, dispers, 100, 950, 450, 525, bottomType = "NONE")
	addComponent(degreepanel, priordegreemean, 500, 475, 950, 50, bottomType = "NONE")
	addComponent(degreepanel, priordegreeSD, 500, 950, 950, 525, bottomType = "NONE")
	#addComponent(dialog, degreedist, 650, 350, 700, 50, bottomType = "NONE")
	
	
	#column 2
	
#prior data panel
	priorpanel <- new(JPanel)
	priorpanel$setBorder(J("javax.swing.BorderFactory")$createTitledBorder("Prior Population Data"))
	priorpanel$setLayout(new(AnchorLayout))
	addComponent(dialog, priorpanel, 375, 950, 900, 500,bottomType="REL")		
	
	addComponent(priorpanel, max_N, 100, 475, 225, 50, bottomType = "REL")
	addComponent(priorpanel, priormed, 100, 950, 225, 525, bottomType = "REL")
	
	addComponent(priorpanel, priormean, 250, 475, 375, 50, bottomType = "REL")
	addComponent(priorpanel, priorsd, 250, 950, 375, 525, bottomType = "REL")
	
	addComponent(priorpanel, priormode, 400, 475, 525, 50, bottomType = "REL")
	addComponent(priorpanel, quarts, 400, 950, 525, 525, bottomType = "REL")
	
	addComponent(priorpanel, priormodeprop, 550, 750, 675, 50, bottomType = "REL")
	#addComponent(priorpanel, priormedprop, 550, 950, 675, 525, bottomType = "REL")	
	addComponent(priorpanel, typedist, 725, 950, 950, 50, bottomType = "REL")

	checkFunc <- function(x) {
		
		if (degreevar$getRModel() == "c()") 
			return("Please enter degree variable")
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
		priorsizedistribution=typedist$getRModel()
		
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
		
		#no entry field
		alpha = "NULL"
		degreedistribution = "\"cmp\"" #
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
		
#		checked to here
		cmd <- "posize <- posteriorsize(" %+% s %+% #", median.prior.size=" %+% median.prior.size %+% ", interval =" %+% intervalsize %+%
					", burnin=" %+% burnin %+%
					", maxN=" %+% maxN %+% ", K=" %+% K %+% 
					", samplesize=" %+% samplesize %+% 
					", quartiles.prior.size=" %+% quartiles.prior.size %+% 
					", mean.prior.size =" %+% mean.prior.size %+% 
					", mode.prior.size=" %+% mode.prior.size %+% 
					", priorsizedistribution=" %+% priorsizedistribution %+% 
					", effective.prior.df= 1" %+% 
					", sd.prior.size =" %+%	sd.prior.size %+% 
					", mode.prior.sample.proportion=" %+% mode.prior.sample.proportion %+% 
					", alpha=" %+% alpha %+% 
					", degreedistribution=" %+% degreedistribution %+% 
					", mean.prior.degree=" %+% mean.prior.degree %+% 
					", sd.prior.degree=" %+% sd.prior.degree %+% 
					", df.mean.prior =" %+% df.mean.prior %+% 
					", df.sd.prior=" %+% df.sd.prior %+% 
					", Np=" %+% Np %+% 
					", dispersion=" %+% dispersion %+% 
					", nk=tabulate(s,nbin=K), n=length(s)" %+% 
					", muproposal=0.1, sigmaproposal=0.15, burnintheta=500" %+% 
					", parallel=1, parallel.type=\"PVM\", seed=NULL, verbose=TRUE" %+% 
					")"
		
			
			#add mean0 and mean1
			#samping defaults are different for posterior size
			
		if(diseasebox$getModel()$size()>0) { #use posteriordisease function				
				
				dis2<- unlist(strsplit(diseasevar$getRModel(), "[\"]"))[2]		
				dis <- varSel$getModel() %+% "$" %+% dis2
				Np0 = "0"
				Np1 = "0"
				burnintheta="500"	
				cmd <- "podisease <- posteriordisease(" %+% s %+% ", dis=" %+% dis %+%
						#", mean0.prior.degree =" %+% mean0.prior.degree %+%
						#", mean1.prior.degree =" %+% mean1.prior.degree %+%
						", sd.prior.degree =" %+% sd.prior.degree %+%
						", df.mean.prior =" %+% df.mean.prior %+% 
						", df.sd.prior=" %+% df.sd.prior %+% 
						", Np0=" %+% 0 %+% 
						", Np1=" %+% 0 %+% 
						", samplesize=" %+% samplesize %+%	
						", burnin=" %+% burnin %+%	
						", interval=" %+% interval %+%	
						", burnintheta=" %+% burnintheta %+%	
						", priorsizedistribution=" %+% priorsizedistribution %+%	
						", mean.prior.size =" %+% mean.prior.size %+% 
						", sd.prior.size =" %+%	sd.prior.size %+% 
						", mode.prior.sample.proportion=" %+% mode.prior.sample.proportion %+% 
						", median.prior.size=" %+% median.prior.size %+% 
						", mode.prior.size =" %+% mode.prior.size %+% 
						", quartiles.prior.size=" %+% quartiles.prior.size %+%
						", effective.prior.df= 1" %+% 
						", alpha=" %+% alpha %+% 
						", degreedistribution=" %+% degreedistribution %+% 
						", maxN=" %+% maxN %+% 
						", K=" %+% K %+% 
						", n=length(s)," %+%
						", dispersion=" %+% dispersion %+%
						", nk0=tabulate(s[dis==0],nbin=K), nk1=tabulate(s[dis==1],nbin=K)" %+%
						", muproposal=0.1, sigmaproposal=0.15, parallel=1, seed=NULL" %+%
						", verbose=TRUE"  %+%
						")"

				cmd <- cmd %+% #";podisease\n"
							}
		execute(cmd)
				
	}	
	dialog$setRunFunction(toJava(runFunc))
	dialog
}