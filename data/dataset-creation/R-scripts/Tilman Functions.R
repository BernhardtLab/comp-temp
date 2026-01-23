#Functions for Tilman script

se <- function(x) sd(x, na.rm=TRUE)/sqrt(length(x[!is.na(x)]))
filtertparameter<-function(d,parameter) filter(d,Grouped.Parameter==parameter)
tdatamean<-function(d,parameter) mean(filtertparameter(d,parameter)$Activation.Energy, na.rm = TRUE)
listtmeans<-function(d) c(tdatamean(d,"k: growth"),
                          tdatamean(d,"k: nutrient uptake"),
                          tdatamean(d,"mortality rate"),
                          tdatamean(d,"growth rate"),
                          tdatamean(d,"consumption rate"))
tdatase<-function(d,parameter) se(filtertparameter(d,parameter)$Activation.Energy)
listtses<-function(d) c(tdatase(d,"k: growth"),
                        tdatase(d,"k: nutrient uptake"),
                        tdatase(d,"mortality rate"),
                        tdatase(d,"growth rate"),
                        tdatase(d,"consumption rate"))
tdataconf<-function(d,parameter) t.test(filtertparameter(d,parameter)$Activation.Energy, na.rm = TRUE)$conf.int[2]-tdatamean(d,parameter)
listtconfs<-function(d) c(tdataconf(d,"k: growth"),
						  tdataconf(d,"k: nutrient uptake"),
						  NA, #no mortality rates, yet, this was the easiest way to avoid errors
						  tdataconf(d,"growth rate"),
						  tdataconf(d,"consumption rate"))
tdatamax<-function(d,parameter) max(filtertparameter(d,parameter)$Activation.Energy)
listtmax<-function(d) c(tdatamax(d,"k: growth"),
                        tdatamax(d,"k: nutrient uptake"),
                        NA,  #no mortality rates, yet, this was the easiest way to avoid -inf results
                        tdatamax(d,"growth rate"),
                        tdatamax(d,"consumption rate"))

create.n<-function(dataframe) n.equals(as.character(dataframe$NumSpecies))
n.equals<-function(string) paste("N",string,sep="=")

#Functions for including Chen et al. (2014) dataset
tgrowthmean<-function(d,chen_version) mean(c(filtertparameter(d,"growth rate")$Activation.Energy,chen_version$Activation.energy..eV.),na.rm = TRUE)
listcmeans<-function(d,chen_version) c(tdatamean(d,"k: growth"),
									   tdatamean(d,"k: nutrient uptake"),
									   tdatamean(d,"mortality rate"),
									   tgrowthmean(d,chen_version),
									   tdatamean(d,"consumption rate"))
tgrowthse<-function(d,chen_version) se(c(filtertparameter(d,"growth rate")$Activation.Energy,chen_version$Activation.energy..eV.))
listcses<-function(d,chen_version) c(tdatase(d,"k: growth"),
									 tdatase(d,"k: nutrient uptake"),
									 tdatase(d,"mortality rate"),
									 tgrowthse(d,chen_version),
									 tdatase(d,"consumption rate"))
tgrowthconf<-function(d,chen_version) t.test(c(chen_version$Activation.energy..eV., filter(d,Grouped.Parameter=="growth rate")$Activation.Energy),na.rm = TRUE)$conf.int[2]-tgrowthmean(d,chen_version)
listcconfs<-function(d,chen_version) c(tdataconf(d,"k: growth"),
									   tdataconf(d,"k: nutrient uptake"),
									   NA, #no mortality rates, yet, this was the easiest way to avoid errors
									   tgrowthconf(d,chen_version),
									   tdataconf(d,"consumption rate"))
tgrowthmax<-function(d,chen_version) max(filtertparameter(d,"growth rate")$Activation.Energy,chen_version$Activation.energy..eV., na.rm = TRUE)
listcmax<-function(d,chen_version) c(tdatamax(d,"k: growth"),
									 tdatamax(d,"k: nutrient uptake"),
									 NA, #no mortality rates, yet, this was the easiest way to avoid -inf results
									 tgrowthmax(d,chen_version),
									 tdatamax(d,"consumption rate"))
tgrowthmin<-function(d,chen_version) min(filtertparameter(d,"growth rate")$Activation.Energy,chen_version$Activation.energy..eV., na.rm = TRUE)
tdatamin<-function(d,parameter) min(filtertparameter(d,parameter)$Activation.Energy)
listcmin<-function(d,chen_version) c(tdatamin(d,"k: growth"),
									 tdatamin(d,"k: nutrient uptake"),
									 NA, #no mortality rates, yet, this was the easiest way to avoid -inf results
									 tgrowthmin(d,chen_version),
									 tdatamin(d,"consumption rate"))