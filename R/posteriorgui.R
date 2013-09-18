# Author: JaneAc

#NOTES/TO DO:

#mention (in tooltip?) that prior data can be calculated via prior distribution dialog?
#add info buttons on subdialogs with screenshots of annotated versions of the dialogs
#nk=tabulate(s,nbin=K)" %+%  caused crash when included excplicitly. also a problem in the marginal functions.
###############################################################################

.makePosteriorDistribution <- function() {
	
	JPanel <- J("javax.swing.JPanel")
	JButton <- J("javax.swing.JButton")
	JSeparator <- J("javax.swing.JSeparator")
	Dimension <- J("java.awt.Dimension")
	AnchorLayout <- J("org.rosuda.JGR.layout.AnchorLayout")
	SingletonDJList <- J("org.rosuda.deducer.toolkit.SingletonDJList")	
	ActionListener <- J("org.rosuda.deducer.widgets.event.RActionListener")	
	
	
#Components
	
	##TOP: choose degree and optional status variables. 
	dialog <- new(Deducer::SimpleRDialog)
	dialog$setSize(500L, 550L)
	dialog$setTitle("Calculate Posterior Distribution of Population Size")
	
	varSel <- new(Deducer::VariableSelectorWidget)
	varSel$setRDataFilter("is.rds.data.frame")
	
	#degreevar <- new(Deducer::SingleVariableWidget,"Degree Variable",varSel)
	diseasevar <- new(Deducer::SingleVariableWidget,"Status (Optional)",varSel)
	wavevar <- new(Deducer::SingleVariableWidget,"Order (Optional)",varSel) #Optional
	
	##choose to use posteriorsize or posteriordisease 
	#use to be a checkbox:
	#diseasebox <- new(Deducer::CheckBoxesWidget,.jarray("Use Status Variable")) 
	#now just use whether status var is filled in
	
	
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
	
	#for posterior disease
	intervalsize2 <- new(Deducer::TextFieldWidget, "Interval")
	intervalsize2$setDefaultModel("1")
	intervalsize2$setLowerBound(1)
	
	burn2 <- new(Deducer::TextFieldWidget, "Burnin")
	burn2$setDefaultModel("100")
	burn2$setLowerBound(0)
	
	samples2 <- new(Deducer::TextFieldWidget, "Number of Samples")
	samples2$setDefaultModel("1000")
	samples2$setLowerBound(1)

	#POP PRIOR:	
	#max_N <- new(Deducer::TextFieldWidget, "Pop. Max")
	#max_N$setLowerBound(1)

	priormed <- new(Deducer::TextFieldWidget, "Median")
	priormed$setLowerBound(1)
	
	priormean <- new(Deducer::TextFieldWidget, "Mean")
	priormean$setLowerBound(1)
	
	priorsd <- new(Deducer::TextFieldWidget, "Standard Deviation (optional)")
	priorsd$setLowerBound(0)
	
	priormode <- new(Deducer::TextFieldWidget, "Mode")
	priormode$setLowerBound(1)
	
	quarts25 <- new(Deducer::TextFieldWidget, "Quartiles - 25%")
	quarts25$setLowerBound(0)
	quarts75 <- new(Deducer::TextFieldWidget, "75%")
	quarts75$setLowerBound(0)
	
	#priormodeprop <- new(Deducer::TextFieldWidget, "Prop. Mode")
	#priormodeprop$setLowerBound(0)
	#priormodeprop$setUpperBound(1)
	#priormodeprop$setDefaultModel(".5")
	
	types = c("Beta","Flat") #"Neg-binom","Poisson-log-norm"
	typedist <-new(Deducer::ComboBoxWidget, types)
	typedist$setTitle("Prior Pop. Size Distribution ", TRUE)
	typedist$setDefaultModel("Beta")
	
	##DEGREE - two copies of some fields because the posteriorsize and posteriordisease functions are run through the same dialog
	maxDeg1 <- new(Deducer::TextFieldWidget, "Max Degree (K)") #= K, default = round(quantile(s,0.80))
	maxDeg1$setLowerBound(1)
	maxDeg2 <- new(Deducer::TextFieldWidget, "Max Degree (K)") #= K, default = round(quantile(s,0.80))
	maxDeg2$setLowerBound(1)
	
	degs1 = c("Conway-Maxwell-Poisson","Neg-binom","Poisson-log-normal")
	degreedist1 <-new(Deducer::ComboBoxWidget, degs1)
	degreedist1$setTitle("Degree Distribution Shape", TRUE)
	degreedist1$setDefaultModel("Conway-Maxwell-Poisson")
	
	degreedist2 <-new(Deducer::ComboBoxWidget, degs1)
	degreedist2$setTitle("Degree Distribution Shape", TRUE)
	degreedist2$setDefaultModel("Conway-Maxwell-Poisson")
	
	priordegreemean <- new(Deducer::TextFieldWidget, "Mean Degree")
	priordegreemean$setLowerBound(0)
	
	priordegreeSD1 <- new(Deducer::TextFieldWidget, "Standard Deviation")
	priordegreeSD1$setLowerBound(0)
	priordegreeSD2 <- new(Deducer::TextFieldWidget, "Standard Deviation")
	priordegreeSD2$setLowerBound(0)
	
	priordegreemean0 <- new(Deducer::TextFieldWidget, "Mean, Status = 0")
	priordegreemean0$setLowerBound(0)
	priordegreemean0$setDefaultModel("7")
	
	priordegreemean1 <- new(Deducer::TextFieldWidget, "Mean, Status = 1")
	priordegreemean1$setLowerBound(0)
	priordegreemean1$setDefaultModel("7")
	
	#dispers <- new(Deducer::TextFieldWidget, "Dispersion")
	#dispers$setLowerBound(0)

	
	#Buttons 
	mcbutton <- new(JButton,"Sampling Prefs")
	mcbutton$setToolTipText("Set burnin period, sampling interval and number of samples for MCMC sampling.")
	
	mcsubdialog <- new(SimpleRSubDialog,dialog,"Set MCMC Sampling Parameters")
	setSize(mcsubdialog,300L,200L)
	#for posteriordisease:
	mcsubdialog2 <- new(SimpleRSubDialog,dialog,"Set MCMC Sampling Parameters")
	setSize(mcsubdialog2,300L,200L)
	
	mcactionFunction <- function(cmd,ActionEvent){
		if (diseasevar$getRModel()=="c()") {mcsubdialog<-mcsubdialog}
			else {mcsubdialog<-mcsubdialog2}	
			mcsubdialog$setLocationRelativeTo(mcbutton)
			mcsubdialog$run()
	}
	mclistener <- new(ActionListener)
	mclistener$setFunction(toJava(mcactionFunction))
	mcbutton$addActionListener(mclistener)
	
	degbutton <- new(JButton,"Prior Degree Knowledge")
	degbutton$setToolTipText("Status = 0 and Status = 1 fields only if \"Use Status Var\" is checked.")
	
	#one sub-dialog for posteriorsize and one for posteriordisease
	degsubdialog1 <- new(SimpleRSubDialog,dialog,"Prior Degree Knowledge")
	degsubdialog2 <- new(SimpleRSubDialog,dialog,"Prior Degree Knowledge")
	
	setSize(degsubdialog1,300L,300L)
	setSize(degsubdialog2,300L,300L)
	
	addComponent(degsubdialog1, priordegreemean, 100, 475, 250, 50, bottomType = "NONE")
	addComponent(degsubdialog1, maxDeg1, 300, 475, 450, 50, bottomType = "NONE")
	addComponent(degsubdialog1, priordegreeSD1, 300, 950, 450, 525, bottomType = "NONE")
	#addComponent(degsubdialog, dispers, 100, 950, 200, 525, bottomType = "NONE")
	addComponent(degsubdialog1, degreedist1, 575, 950, 800, 50, bottomType = "REL")
	
	addComponent(degsubdialog2, priordegreemean0, 100, 475, 250, 50, bottomType = "NONE")
	addComponent(degsubdialog2, priordegreemean1, 100, 950, 250, 525, bottomType = "NONE")
	addComponent(degsubdialog2, maxDeg2, 300, 475, 450, 50, bottomType = "NONE")
	addComponent(degsubdialog2, priordegreeSD2, 300, 950, 450, 525, bottomType = "NONE")
	addComponent(degsubdialog2, degreedist2, 575, 950, 800, 50, bottomType = "REL")

	
	degactionFunction <- function(cmd,ActionEvent){
		if (diseasevar$getRModel()=="c()") {degsubdialog<-degsubdialog1}
		else {degsubdialog<-degsubdialog2}	
		degsubdialog$setLocationRelativeTo(degbutton)
		degsubdialog$run()
	}
	deglistener <- new(ActionListener)
	deglistener$setFunction(toJava(degactionFunction))
	degbutton$addActionListener(deglistener)
	
	popbutton <- new(JButton,"Prior Pop. Size Knowledge")
	
	popsubdialog <- new(SimpleRSubDialog,dialog,"Prior Population Knowledge")
	setSize(popsubdialog,300L,500L)
	
	popactionFunction <- function(cmd,ActionEvent){
		popsubdialog$setLocationRelativeTo(popbutton)
		popsubdialog$run()
	}
	poplistener <- new(ActionListener)
	poplistener$setFunction(toJava(popactionFunction))
	popbutton$addActionListener(poplistener)
	
	plotbox <- new(Deducer::CheckBoxesWidget,.jarray("Plot Distribution"))
	setSize(plotbox,400,75)	
	plotbox$setDefaultModel(c("Plot Distribution"))
	
	
	
#Arrange Components:

	addComponent(dialog, varSel, 50, 450, 500, 50)
	
	addComponent(dialog, plotbox, 60, 950, 150, 620)
	addComponent(dialog, diseasevar, 200, 950, 350, 500, topType="REL",bottomType="REL")
	addComponent(dialog, wavevar, 375, 950, 525, 500, topType="REL",bottomType="REL")
	

	#buttons
	addComponent(dialog,mcbutton,590,450,690,50)
	addComponent(dialog,degbutton,590,950,690,550)
	addComponent(dialog, typedist, 750,450,850,50, bottomType = "REL")
	addComponent(dialog,popbutton,750,950,850,550)	
	
	addComponent(mcsubdialog, burn, 100, 475, 450, 50, bottomType = "NONE")
	addComponent(mcsubdialog, intervalsize, 100, 950, 450, 525, bottomType = "NONE")
	addComponent(mcsubdialog, samples, 500, 950, 850, 50, bottomType = "NONE")

	addComponent(mcsubdialog2, burn2, 100, 475, 450, 50, bottomType = "NONE")
	addComponent(mcsubdialog2, intervalsize2, 100, 950, 450, 525, bottomType = "NONE")
	addComponent(mcsubdialog2, samples2, 500, 950, 850, 50, bottomType = "NONE")
	
	#column 2
	
	mmmpanel <- new(JPanel)
	mmmpanel$setBorder(J("javax.swing.BorderFactory")$createTitledBorder("Fill in one row or none (default)"))
	mmmpanel$setLayout(new(AnchorLayout))
	addComponent(popsubdialog, mmmpanel, 50, 950, 700, 50,bottomType="REL")		

	#degreeinfo
	
	addComponent(mmmpanel, priormean, 100, 950, 250, 50, bottomType = "NONE")
	addComponent(mmmpanel, priormed, 300, 950, 450, 50, bottomType = "NONE")
	addComponent(mmmpanel, priormode, 500, 950, 650, 50, bottomType = "NONE")
	addComponent(mmmpanel, quarts25, 700, 490, 850, 50, bottomType = "NONE")
	addComponent(mmmpanel, quarts75, 700, 950, 850, 510, bottomType = "NONE")
	addComponent(popsubdialog, priorsd, 725, 900, 840, 100, bottomType = "REL")	

#Functions
	
	checkFunc <- function(x) {
		
		if (varSel$getRModel()=="c()") {
			return("Please select a data set")}
		if (xor(quarts25$getModel()=="",quarts75$getModel()==""))
			return("Quartile entries should be empty or both filled")
			
		else("")
		
	}
	dialog$setCheckFunction(toJava(checkFunc))	
	
	
	#if(plotbox$getModel()$size()>0)
	runFunc <- function(x){
		
		dataset <- unlist(strsplit(varSel$getRModel(), "[\"]"))[2]
		"%+%" <- function(x, y) paste(x, y, sep = "")
		
		if (wavevar$getRModel() != "c()") {
			ordername <- unlist(strsplit(wavevar$getRModel(), "[\"]"))[2]
			order.var <- dataset %+% "$" %+% ordername
			s <- "get.net.size(" %+% dataset %+% ")[order(" %+% order.var %+% ")]"		
		}		
		
		else {s <- "get.net.size(" %+% dataset %+% ")[order(get.wave(" %+% dataset %+% "))]"}#result from get.wave shouldn't have na's
		
		priorsizedistribution <- switch(typedist$getModel(),
				"Beta"="beta",
				#"Neg-binom"="nbinom",
				#"Poisson-log-norm"="pln",
				"Flat" = "flat" )
		
		quarts1 <- quarts25$getRModel()
		quarts2 <- quarts75$getRModel()
		
		cmd2 <- ""
		#pop fields
		if (priormed$getModel()!="") {
			med.prior.size <- priormed$getModel()
			cmd2 <- cmd2  %+% ", med.prior.size =" %+% med.prior.size
		}	
		
		if (priormean$getModel()!="") {
			mean.prior.size <- priormean$getModel()
			cmd2 <- cmd2  %+% ", mean.prior.size =" %+% mean.prior.size
		}
		
		if (priormode$getModel()!="") {
			mode.prior.size <- priormode$getModel()
			cmd2 <- cmd2  %+% ", mode.prior.size =" %+% mode.prior.size
		}
		
		if (priorsd$getModel()!="") {
			sd.prior.size = priorsd$getModel()
			cmd2 <- cmd2  %+% ", sd.prior.size =" %+% sd.prior.size
		}
		
		if (quarts1!="\"\"" && quarts2!="\"\"") {#doublechecks that both must be filled in if one is filled in
			quartiles.prior.size = c(quarts1,quarts2)
			cmd2 <- cmd2 %+% ", quartiles.prior.size=" %+% quartiles.prior.size
		}
		
	
#### OPTION 1		
		
		if (diseasevar$getRModel()=="c()") {
			
		priordegreedistribution <- switch(degreedist1$getModel(),
				"Conway-Maxwell-Poisson" = "cmp", 
				"Neg-binom"="nbinom",
				"Poisson-log-norm"="pln")
			
		cmd <- "posize <- posteriorsize(s =" %+% s %+% 
				   #", median.prior.size=" %+% median.prior.size %+% ", 
				   #", maxN=" %+% maxN %+%
					", priorsizedistribution= \"" %+% priorsizedistribution %+% "\""
		
	#add to command depending on input:
		
			#degreefields
			if (priordegreemean$getModel()!="") {
				mean.prior.degree = priordegreemean$getModel()
				cmd <- cmd  %+% ", mean.prior.degree=" %+% mean.prior.degree		
			}
		
			if (maxDeg1$getModel()!="") {
				K = maxDeg1$getModel()
				cmd <- cmd  %+% ", K =" %+% K
			}
		
			if (priordegreeSD1$getModel()!="") {
				sd.prior.degree = priordegreeSD1$getModel()
				cmd <- cmd  %+% ", sd.prior.degree=" %+% sd.prior.degree 
			}
		
			if (samples$getModel()!="" && samples$getModel()!="1000") {
				samplesize <- samples$getModel()
				cmd <- cmd  %+% ", samplesize= " %+% samplesize
			}
			
			if (burn$getModel()!="" && burn$getModel()!="5000") {
				burnin <- burn$getModel()
				cmd <- cmd  %+% ", burnin= " %+% burnin
			}
			
			if (intervalsize$getModel()!="" && intervalsize$getModel()!="10") {
				interval <- intervalsize$getModel()
				cmd <- cmd  %+% ", interval = " %+% interval
			}

		cmd <-cmd %+% cmd2 # %+% "\n posize\n" #print this?

		cmd <-cmd %+% ")"

		if(plotbox$getModel()$size()>0) {cmd <- cmd %+% "\nplot.size(posize)"}
		
		}

#### OPTION 2		
				
		else { #use posteriordisease function				
				
			priordegreedistribution <- switch(degreedist2$getModel(),
					"Conway-Maxwell-Poisson" = "cmp", 
					"Neg-binom"="nbinom",
					"Poisson-log-norm"="pln")
			
			if (wavevar$getRModel() != "c()") {
				ordername <- unlist(strsplit(wavevar$getRModel(), "[\"]"))[2]
				order.var <- dataset %+% "$" %+% ordername
				s <- "get.net.size(" %+% dataset %+% ")[order(" %+% order.var %+% ")]"		
			}		
			
			else {s <- "get.net.size(" %+% dataset %+% ")[order(get.wave(" %+% dataset %+% "))]"}#result from get.wave shouldn't have na's
			
			
			dis2<- unlist(strsplit(diseasevar$getRModel(), "[\"]"))[2]		
			dis <- dataset %+% "$" %+% dis2
			
			if (wavevar$getRModel() != "c()") {
				ordername <- unlist(strsplit(wavevar$getRModel(), "[\"]"))[2]
				order.var <- dataset %+% "$" %+% ordername
				dis <- dis %+% "[order(" %+% order.var %+% ")]"
			}		
			
			else {disease.var <- dis %+% "[order(get.wave(" %+% dataset %+% "))]"} #result from get.wave shouldn't have na's

			mean0.prior.degree = priordegreemean0$getModel()
			mean1.prior.degree = priordegreemean1$getModel()			
			
			cmd <- "podisease <- posteriordisease(s = " %+% s %+% 
						", dis = " %+% dis

				if (mean0.prior.degree!="" && mean0.prior.degree!="7") {
					mean0.prior.degree <- priordegreemean0$getModel()
					cmd <- cmd  %+% ", mean0.prior.degree = " %+% mean0.prior.degree
					}
					
				if (mean1.prior.degree!="" && mean1.prior.degree!="7") {
					mean1.prior.degree <- priordegreemean0$getModel()
					cmd <- cmd  %+% ", mean1.prior.degree = " %+% mean1.prior.degree
					}
					
				if (priordegreeSD2$getModel()!="") {
					sd.prior.degree = priordegreeSD2$getModel()
					cmd <- cmd  %+% ", sd.prior.degree=" %+% sd.prior.degree 
					}
					
				if (maxDeg2$getModel()!="") {
					K = maxDeg2$getModel()
					cmd <- cmd  %+% ", K =" %+% K
					}
					
				if (samples2$getModel()!="" && samples2$getModel()!="1000") {
					samplesize <- samples2$getModel()
					cmd <- cmd  %+% ", samplesize= " %+% samplesize
					}
					
				if (burn2$getModel()!="" && burn2$getModel()!="100") {
					burnin <- burn2$getModel()
					cmd <- cmd  %+% ", burnin= " %+% burnin
					}
					
				if (intervalsize2$getModel()!="" && intervalsize2$getModel()!="1") {
					interval <- intervalsize2$getModel()
					cmd <- cmd  %+% ", interval = " %+% interval
					}
					
				cmd <- cmd %+% ", degreedistribution= \"" %+% priordegreedistribution %+% "\""  
				cmd <- cmd %+% cmd2 %+% ")"

				if(plotbox$getModel()$size()>0) {cmd <- cmd %+% "\nplot.size(podisease)"}
				
						}
		execute(cmd)

	}	
	dialog$setRunFunction(toJava(runFunc))
	dialog
}
