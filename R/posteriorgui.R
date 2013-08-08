# TODO: Add comment
# 
# Author: JaneAc
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
	#choose to use posteriorsize or posterior disease function via checkbox (diseasebox)
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
	
	burnin <- new(Deducer::TextFieldWidget, "Burnin")
	burnin$setDefaultModel("5000")
	burnin$setLowerBound(0)
	
	samples <- new(Deducer::TextFieldWidget, "Number of samples")
	samples$setDefaultModel("1000")
	samples$setLowerBound(1)
	
	
	#column 2
	#prior information:
	
	#can be calculated via prior distribution dialog
	
	maxN <- new(Deducer::TextFieldWidget, "Population Max")
	maxN$setLowerBound(1)

	priormed <- new(Deducer::TextFieldWidget, "Prior Median")
	priormed$setLowerBound(1)
	
	priormean <- new(Deducer::TextFieldWidget, "Prior Mean")
	priormean$setLowerBound(1)
	
	priorsd <- new(Deducer::TextFieldWidget, "Prior SD")
	priorsd$setLowerBound(0)
	
	priormode <- new(Deducer::TextFieldWidget, "Prior Mode")
	priormode$setLowerBound(1)
	
	quarts <- new(Deducer::TextFieldWidget, "Quartiles (25%,75%)")
	
	priormodeprop <- new(Deducer::TextFieldWidget, "Proportion Prior Mode")
	priormodeprop$setLowerBound(0)
	priormodeprop$setUpperBound(1)
	priormodeprop$setDefaultModel(".5")
	
	#MAKE THIS A DROP DOWN
	typedist <-new(Deducer::ButtonGroupWidget, "Prior Dist Type", c("proportion","flat","nbinom","pln"))
	typedist$setDefaultModel("proportion")
	
	#degree info
	
	maxDeg <- new(Deducer::TextFieldWidget, "Max Degree") #= K, default = round(quantile(s,0.80))
	maxDeg$setLowerBound(1)
	maxDeg$setToolTipText("The maximum degree for an individual. This is usually calculated as twice the maximum observed degree.")
	
	degreedist <-new(Deducer::ButtonGroupWidget, "Degree Distribution Type", c("cmp","nbinom","pln"))
	degreedist$setDefaultModel("CMP")
	#below isn't working. Doesn't work on ButtonGroupWidget?
	#degreedist$setToolTipText("cmp = Conway-Maxwell-Poisson; nbinom = Negative Binomial; pln = Poisson-log-normal")
	
	priordegreemean <- new(Deducer::TextFieldWidget, "Degree Mean")
	priordegreemean$setLowerBound(0)
	
	priordegreeSD <- new(Deducer::TextFieldWidget, "Degree SD")
	priordegreeSD$setLowerBound(0)
	
	dispersion <- new(Deducer::TextFieldWidget, "Dispersion")
	dispersion$setLowerBound(0)
	
	
	#plotbox <- new(Deducer::CheckBoxesWidget,.jarray("Plot Distribution"))
	#setSize(plotbox,400,75)	
	#plotbox$setDefaultModel(c("Plot Distribution"))
	
	#column2
	priormean <- new(Deducer::TextFieldWidget, "Prior Mean")
	priormean$setLowerBound(1)

	
	priormode <- new(Deducer::TextFieldWidget, "Prior Mode")
	priormed$setLowerBound(1)
	
	priorquartiles <- new(Deducer::TextFieldWidget, "Prior Quartiles (25%, 75%)")
	#set to require acoordinate pair
	
	#ignore for now - user can set on command line
	#prioralpha <- new(Deducer::TextFieldWidget, "Alpha")
	#effectivedf <- new(Deducer::TextFieldWidget, "Effective Prior DF")
	#effectivedf$setDefaultModel("1")
	
	priormedprop <- new(Deducer::TextFieldWidget, "Proportion Prior Median")
	priormedprop$setLowerBound(0)
	priormedprop$setUpperBound(1)
	
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
	
	addComponent(samplepanel, burnin, 100, 475, 450, 50, bottomType = "NONE")
	addComponent(samplepanel, intervalsize, 100, 950, 450, 525, bottomType = "NONE")
	addComponent(samplepanel, samples, 500, 950, 950, 50, bottomType = "NONE")
	
	#degreeinfo
	degreepanel <- new(JPanel)
	degreepanel$setBorder(J("javax.swing.BorderFactory")$createTitledBorder("Prior Degree Data"))
	degreepanel$setLayout(new(AnchorLayout))
	addComponent(dialog, degreepanel, 650, 450, 900, 50,bottomType="REL")		

	addComponent(degreepanel, maxDeg, 100, 475, 450, 50, bottomType = "NONE")
	addComponent(degreepanel, dispersion, 100, 950, 450, 525, bottomType = "NONE")
	addComponent(degreepanel, priordegreemean, 500, 475, 950, 50, bottomType = "NONE")
	addComponent(degreepanel, priordegreeSD, 500, 950, 950, 525, bottomType = "NONE")
	#addComponent(dialog, degreedist, 650, 350, 700, 50, bottomType = "NONE")
	
	
	#column 2
	
#prior data panel
	priorpanel <- new(JPanel)
	priorpanel$setBorder(J("javax.swing.BorderFactory")$createTitledBorder("Prior Pop. Data"))
	priorpanel$setLayout(new(AnchorLayout))
	addComponent(dialog, priorpanel, 375, 950, 900, 500,bottomType="REL")		
	
	addComponent(priorpanel, maxN, 100, 475, 225, 50, bottomType = "REL")
	addComponent(priorpanel, priormed, 100, 950, 225, 525, bottomType = "REL")
	
	addComponent(priorpanel, priormean, 250, 475, 375, 50, bottomType = "REL")
	addComponent(priorpanel, priorsd, 250, 950, 375, 525, bottomType = "REL")
	
	addComponent(priorpanel, priormode, 400, 475, 525, 50, bottomType = "REL")
	addComponent(priorpanel, quarts, 400, 950, 525, 525, bottomType = "REL")
	
	addComponent(priorpanel, priormodeprop, 550, 475, 675, 50, bottomType = "NONE")
	addComponent(priorpanel, priormedprop, 550, 950, 675, 525, bottomType = "NONE")
	
	#Left space to add this once I make it a drop-down. Also, dist for degree?
	#addComponent(dialog, typedist, 900, 900, 1000, 600, bottomType = "NONE")

	checkFunc <- function(x) {}
	dialog$setCheckFunction(toJava(checkFunc))	
	
	runFunc <- function(x){}	
	dialog$setRunFunction(toJava(runFunc))
	dialog
}