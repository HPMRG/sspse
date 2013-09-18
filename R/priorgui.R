# TODO: Add comment
# 
# Author: JaneAc
###############################################################################


.makePriorDistribution <- function() {

	JPanel <- J("javax.swing.JPanel")
	Dimension <- J("java.awt.Dimension")
	AnchorLayout <- J("org.rosuda.JGR.layout.AnchorLayout")
	ActionListener <- J("org.rosuda.deducer.widgets.event.RActionListener")	
	
	#column1	
	
	dialog <- new(Deducer::SimpleRDialog)
	dialog$setSize(480L, 420L)
	dialog$setTitle("Calculate Prior Distribution of Population Size")
	
	samplesize <- new(Deducer::TextFieldWidget, "Sample Size")
	samplesize$setInteger(TRUE)
	samplesize$setLowerBound(1)
	
	maxN <- new(Deducer::TextFieldWidget, "Population Max")
	maxN$setInteger(TRUE)
	maxN$setLowerBound(1)
	
	typedist <-new(Deducer::ButtonGroupWidget, "Distribution Type", c("Beta","Flat")) #"Neg-binom","Continuous","Poisson-log-norm",
	#typedist$setSize(400,250)
	
	plotbox <- new(Deducer::CheckBoxesWidget,.jarray("Plot Distribution"))
	setSize(plotbox,400,75)	
	plotbox$setDefaultModel(c("Plot Distribution"))
	
	#column2
	priormean <- new(Deducer::TextFieldWidget, "Mean")
	priormean$setLowerBound(1)
	
	priorsd <- new(Deducer::TextFieldWidget, "Standard Deviation (optional)")
	priorsd$setLowerBound(0)
	
	priormed <- new(Deducer::TextFieldWidget, "Median")
	priormed$setLowerBound(1)
	
	priormode <- new(Deducer::TextFieldWidget, "Mode")
	priormed$setLowerBound(1)

	quarts25 <- new(Deducer::TextFieldWidget, "Quartiles 25%")
	quarts25$setLowerBound(0)
	quarts75 <- new(Deducer::TextFieldWidget, "75%")
	quarts75$setLowerBound(0)
	
	#ignore for now - user can set inconsole
	#prioralpha <- new(Deducer::TextFieldWidget, "Alpha")
	#effectivedf <- new(Deducer::TextFieldWidget, "Effective Prior DF")
	#effectivedf$setDefaultModel("1")
	
	#priormodeprop <- new(Deducer::TextFieldWidget, "Proportion Mode")
	#priormodeprop$setLowerBound(0)
	#priormodeprop$setUpperBound(1)
	#priormodeprop$setDefaultModel(".5")
	
	#priormedprop <- new(Deducer::TextFieldWidget, "Proportion Median")
	#priormedprop$setLowerBound(0)
	#priormedprop$setUpperBound(1)
	
	#left column 1
	addComponent(dialog, samplesize, 50, 450, 200, 100, bottomType = "NONE")
    addComponent(dialog, maxN, 225, 450, 375, 100, bottomType = "NONE")
	addComponent(dialog, typedist, 410, 450, 700, 100, bottomType = "REL")
	addComponent(dialog, plotbox, 700, 450, 900, 100, bottomType = "NONE")

	#right column 2
	
	mmmpanel <- new(JPanel)
	mmmpanel$setBorder(J("javax.swing.BorderFactory")$createTitledBorder("Fill at most one row"))
	mmmpanel$setLayout(new(AnchorLayout))
	
	addComponent(dialog, mmmpanel, 50, 950, 700, 500,bottomType="REL")		
	addComponent(mmmpanel, priormean, 100, 900, 250, 50, bottomType = "NONE")
	addComponent(mmmpanel, priormed, 300, 900, 450, 50, bottomType = "NONE")
	addComponent(mmmpanel, priormode, 500, 900, 650, 50,  bottomType = "NONE")
	addComponent(mmmpanel, quarts25, 700, 455, 850, 50, bottomType = "NONE")
	addComponent(mmmpanel, quarts75, 700, 900, 850, 465, bottomType = "NONE")
	addComponent(dialog, priorsd, 725, 910, 875, 520, bottomType = "NONE")
	
	
	#addComponent(dialog, priormodeprop, 675, 900, 775, 500, bottomType = "NONE")
	#addComponent(dialog, priormedprop, 800, 900, 900, 500, bottomType = "NONE")
	
	checkFunc <- function(x) {
		if (samplesize$getModel() == "") 
			return("Please enter sample size")
		if (maxN$getModel() == "" && typedist$getModel()!="Proportion") 
			return("Please enter population max")
		#if (priorquartiles$getModel()!="" && !length(strsplit(priorquartiles$getModel(),",")[[1]])%in%c(0,2))
		#	return("Quartile entry should be empty or of form low,high e.g 2000,5000")
			#return(strsplit(priorquartiles$getModel(),",")[[1]])
		else("")
}

	dialog$setCheckFunction(toJava(checkFunc))	

	runFunc <- function(x){
	    
		"%+%" <- function(x, y) paste(x, y, sep = "")
		
		n <- samplesize$getModel()
		max_N <- maxN$getModel()
		dist_type <- switch(typedist$getModel(),
				"Beta"="beta",
				#"Neg-binom"="nbinom",
				#"Continuous"="continuous",
				#"Poisson-log-norm"="pln",
				"Flat" = "flat" )
		
		quarts1 <- quarts25$getRModel()
		quarts2 <- quarts75$getRModel()
			
		cmd <- "dsp <- dsizeprior(" %+% n %+%
				", maxN=" %+% max_N %+%
				", type=\"" %+% dist_type %+% "\""
		
		if (priormean$getModel()!="") {
			mean.prior.size = priormean$getModel()
			cmd <- cmd  %+% ", mean.prior.size=" %+% mean.prior.size
		}
		
		if (priormed$getModel()!="") {
			median.prior.size = priormed$getModel()
			cmd <- cmd  %+% ", median.prior.size=" %+% median.prior.size 
		}
				
		if (priormode$getModel()!="") {
			mode.prior.size = priormode$getModel()
			cmd <- cmd  %+% ", mode.prior.size=" %+% mode.prior.size 
		}
				
		if (quarts1!="\"\"" && quarts2!="\"\"") {#doublechecks that both must be filled in if one is filled in
			quartiles.prior.size = "c(" %+% strsplit(quarts1,"\"")[[1]][2] %+% "," %+% strsplit(quarts1,"\"")[[1]][2] %+% ")"
			print(quartiles.prior.size)
			cmd <- cmd %+% ", quartiles.prior.size=" %+% quartiles.prior.size
		}
			
		if (priorsd$getModel()!="") {
			sd.prior.size = priorsd$getModel()
			cmd <- cmd  %+% ", sd.prior.size=" %+% sd.prior.size
		}
		
		
		
	if(plotbox$getModel()$size()>0) { #return info and plot instead of pmf vectors
			cmd <- cmd %+% ")\n dsp[3:14] \n plot(dsp$x,dsp$lprior,main = \"Prior Distribution for Population Size (" %+% dist_type %+% ")\", xlab = \"population size\", ylab = \"prior density\", type=\"n\")"
			cmd <- cmd %+% "\nlines(dsp$x,dsp$lprior,lty = 2)"
			}
			
	else (cmd <- cmd %+% ")\n dsp\n")
	
	print(cmd)
	execute(cmd)
	}	
	dialog$setRunFunction(toJava(runFunc))
	dialog
}
