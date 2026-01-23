#Functions for MacArthur script

#Basic Functions
se <- function(x) sd(x, na.rm=TRUE)/sqrt(length(x[!is.na(x)]))
mdatamean<-function(d,parameter) mean(filterparameter(d,parameter)$Activation.energy, na.rm = TRUE)
filterparameter<-function(d,parameter) filter(d,Simple.Parameter==parameter)
listmmeans<-function(d) c(mdatamean(d,"resource growth rate"),
                          mdatamean(d,"resource carrying capacity"),
                          mdatamean(d,"mortality rate"),
                          mdatamean(d,"conversion efficiency"),
                          mdatamean(d,"consumption rate"))
mdatase<-function(d,parameter) se(filterparameter(d,parameter)$Activation.energy)
listmses<-function(d) c(mdatase(d,"resource growth rate"),
                        mdatase(d,"resource carrying capacity"),
                        mdatase(d,"mortality rate"),
                        mdatase(d,"conversion efficiency"),
                        mdatase(d,"consumption rate"))
mdataconf<-function(d,parameter) t.test(filterparameter(d,parameter)$Activation.energy, na.rm = TRUE)$conf.int[2]-mdatamean(d,parameter)
listmconfs<-function(d) c(mdataconf(d,"resource growth rate"),
						  mdataconf(d,"resource carrying capacity"),
						  mdataconf(d,"mortality rate"),
						  mdataconf(d,"conversion efficiency"),
						  mdataconf(d,"consumption rate"))
mdatamax<-function(d,parameter) max(filterparameter(d,parameter)$Activation.energy)
listmmax<-function(d) c(mdatamax(d,"resource growth rate"),
                        mdatamax(d,"resource carrying capacity"),
                        mdatamax(d,"mortality rate"),
                        mdatamax(d,"conversion efficiency"),
                        mdatamax(d,"consumption rate"))

create.n<-function(dataframe) n.equals(as.character(dataframe$NumSpecies))
n.equals<-function(string) paste("N",string,sep="=")

#Functions to separate autotrophs and heterotrophs
ahdatamean<-function(d,parameter,trophy) mean(filterahparameter(d,parameter,trophy)$Activation.energy, na.rm = TRUE)
filterahparameter<-function(d,parameter,trophy) filter(d,Simple.Parameter==parameter&Heterotroph.or.Autotroph==trophy)
listahmeans<-function(d) c(ahdatamean(d,"resource growth rate","heterotroph"),
                           ahdatamean(d,"resource growth rate","autotroph"),
                           ahdatamean(d,"resource carrying capacity","heterotroph"),
                           ahdatamean(d,"resource carrying capacity","autotroph"),
                           ahdatamean(d,"mortality rate","heterotroph"),
                           ahdatamean(d,"conversion efficiency","heterotroph"),
                           ahdatamean(d,"consumption rate","heterotroph"))
ahdatase<-function(d,parameter,trophy) se(filterahparameter(d,parameter,trophy)$Activation.energy)
listahses<-function(d) c(ahdatase(d,"resource growth rate","heterotroph"),
                         ahdatase(d,"resource growth rate","autotroph"),
                         ahdatase(d,"resource carrying capacity","heterotroph"),
                         ahdatase(d,"resource carrying capacity","autotroph"),
                         ahdatase(d,"mortality rate","heterotroph"),
                         ahdatase(d,"conversion efficiency","heterotroph"),
                         ahdatase(d,"consumption rate","heterotroph"))
ahdataconf<-function(d,parameter,trophy) t.test(filterahparameter(d,parameter,trophy)$Activation.energy, na.rm = TRUE)$conf.int[2]-ahdatamean(d,parameter,trophy)
listahconfs<-function(d) c(ahdataconf(d,"resource growth rate","heterotroph"),
						   ahdataconf(d,"resource growth rate","autotroph"),
						   ahdataconf(d,"resource carrying capacity","heterotroph"),
						   ahdataconf(d,"resource carrying capacity","autotroph"),
						   ahdataconf(d,"mortality rate","heterotroph"),
						   ahdataconf(d,"conversion efficiency","heterotroph"),
						   ahdataconf(d,"consumption rate","heterotroph"))
ahdatamax<-function(d,parameter,trophy) max(filterahparameter(d,parameter,trophy)$Activation.energy)
listahmax<-function(d) c(ahdatamax(d,"resource growth rate","heterotroph"),
                         ahdatamax(d,"resource growth rate","autotroph"),
                         ahdatamax(d,"resource carrying capacity","heterotroph"),
                         ahdatamax(d,"resource carrying capacity","autotroph"),
                         ahdatamax(d,"mortality rate","heterotroph"),
                         ahdatamax(d,"conversion efficiency","heterotroph"),
                         ahdatamax(d,"consumption rate","heterotroph"))
ahdatamin<-function(d,parameter,trophy) min(filterahparameter(d,parameter,trophy)$Activation.energy)
listahmin<-function(d) c(ahdatamin(d,"resource growth rate","heterotroph"),
                         ahdatamin(d,"resource growth rate","autotroph"),
                         ahdatamin(d,"resource carrying capacity","heterotroph"),
                         ahdatamin(d,"resource carrying capacity","autotroph"),
                         ahdatamin(d,"mortality rate","heterotroph"),
                         ahdatamin(d,"conversion efficiency","heterotroph"),
                         ahdatamin(d,"consumption rate","heterotroph"))

#Functions for multiple resource plot
count.c.observations<-function(species) length(filter(SingleSpeciesM,Species==species&Simple.Parameter=="consumption rate")$Species)
count.resource.species<-function(consumer_species) length(unique(filter(SingleSpeciesM,Species==consumer_species&Simple.Parameter=="consumption rate")$Resource.Species))

list.o.counts<-function(data) {observation_counts=c()
for (val in data$Species) {
  observation_counts=append(observation_counts,count.c.observations(val))}
observation_counts}

list.r.counts<-function(data) {resource_counts=c()
for (val in data$Species) {
  resource_counts=append(resource_counts, count.resource.species(val))}
resource_counts}

#Functions for including Chen et al. (2014) dataset
ahgrowthmean<-function(d,chen_version) mean(c(chen_version$Activation.energy..eV., filter(d,Simple.Parameter=="resource growth rate"&Heterotroph.or.Autotroph=="autotroph")$Activation.energy),na.rm = TRUE)
listahcmeans<-function(d,chen_version) c(ahdatamean(d,"resource growth rate","heterotroph"),
										 ahgrowthmean(d,chen_version),
										 ahdatamean(d,"resource carrying capacity","heterotroph"),
										 ahdatamean(d,"resource carrying capacity","autotroph"),
										 ahdatamean(d,"mortality rate","heterotroph"),
										 ahdatamean(d,"conversion efficiency","heterotroph"),
										 ahdatamean(d,"consumption rate","heterotroph"))
ahgrowthconf<-function(d,chen_version) t.test(c(chen_version$Activation.energy..eV., filter(d,Simple.Parameter=="resource growth rate"&Heterotroph.or.Autotroph=="autotroph")$Activation.energy),na.rm = TRUE)$conf.int[2]-ahgrowthmean(d,chen_version)
listahcconfs<-function(d,chen_version) c(ahdataconf(d,"resource growth rate","heterotroph"),
										 ahgrowthconf(d,chen_version),
										 ahdataconf(d,"resource carrying capacity","heterotroph"),
										 ahdataconf(d,"resource carrying capacity","autotroph"),
										 ahdataconf(d,"mortality rate","heterotroph"),
										 ahdataconf(d,"conversion efficiency","heterotroph"),
										 ahdataconf(d,"consumption rate","heterotroph"))
