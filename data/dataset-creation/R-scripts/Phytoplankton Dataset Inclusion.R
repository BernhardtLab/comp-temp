#Script for working with Chen et al. (2014) data
library(ggplot2)
library(ggpubr)
library(dplyr)
library(viridis)

#MacArthur Model
source("R-scripts/MacArthur Functions.R")

MacArthur<-read.csv("data_raw/macarthur_v17.csv")
ChenData<-read.csv("data_raw/Chen_dataset_v1.csv") #v2 is unecessary unless there is a need to combine species counts, better to use v1
ChenData$Simple.Parameter<-rep("resource growth rate", length(ChenData$Activation.energy..eV.))

#Subset of the data containing E estimates for individual species only
#(from lab experiments with >2 temperature treatments)
SingleSpeciesM<-filter(MacArthur,Single.species.or.multiple=="single species"&Experiment.Type=="Lab"&Number.of.Temps.Measured>2)
ChenSubset<-filter(ChenData,is.na(Standard.error.of.activation.energy)==FALSE) #according to the text, activation energies where SE is NA are calculated from only two points

#Count Species
NStrains<-length(unique(ChenSubset$Species))
NHeterotrophGrowth<-length(unique(filter(SingleSpeciesM,Simple.Parameter=="resource growth rate"&Heterotroph.or.Autotroph=="heterotroph")$Species))
NAutotrophGrowth<-length(unique(filter(SingleSpeciesM,Simple.Parameter=="resource growth rate"&Heterotroph.or.Autotroph=="autotroph")$Species))
NHeterotrophK<-length(unique(filter(SingleSpeciesM,Simple.Parameter=="resource carrying capacity"&Heterotroph.or.Autotroph=="heterotroph")$Species))
NAutotrophK<-length(unique(filter(SingleSpeciesM,Simple.Parameter=="resource carrying capacity"&Heterotroph.or.Autotroph=="autotroph")$Species))
NMortality<-length(unique(filter(SingleSpeciesM,Simple.Parameter=="mortality rate")$Species))
NConsumption<-length(unique(filter(SingleSpeciesM,Simple.Parameter=="consumption rate")$Species))
NEfficiency<-length(unique(filter(SingleSpeciesM,Simple.Parameter=="conversion efficiency")$Species))

AHSingleSpeciesCMeans<-data.frame("Parameter"=c("resource growth rate",
												"resource growth rate",
												"resource carrying capacity",
												"resource carrying capacity",
												"mortality rate",
												"conversion efficiency",
												"consumption rate"))
AHSingleSpeciesCMeans$Heterotroph.or.Autotroph<-c("heterotroph",
												  "autotroph",
												  "heterotroph",
												  "autotroph",
												  "heterotroph",
												  "heterotroph",
												  "heterotroph")
AHSingleSpeciesCMeans$Mean<-listahcmeans(SingleSpeciesM,ChenSubset)
AHSingleSpeciesCMeans$ConfidenceInterval<-listahcconfs(SingleSpeciesM,ChenSubset)
AHSingleSpeciesCMeans$NumSpecies<-c(NHeterotrophGrowth,NAutotrophGrowth,NHeterotrophK,NAutotrophK,NMortality,NEfficiency,NConsumption)
AHSingleSpeciesCMeans$N<-create.n(AHSingleSpeciesCMeans)
AHSingleSpeciesCMeans$Max<-listahmax(SingleSpeciesM)
AHSingleSpeciesCMeans$Min<-listahmin(SingleSpeciesM)
AHSingleSpeciesCMeans$Height<-c(AHSingleSpeciesCMeans$Max[1]+0.15,
								AHSingleSpeciesCMeans$Min[2]-0.15,
								AHSingleSpeciesCMeans$Min[3]-0.15,
								AHSingleSpeciesCMeans$Max[4]+0.15,
								AHSingleSpeciesCMeans$Max[5]+0.15,
								AHSingleSpeciesCMeans$Max[6]+0.15,
								AHSingleSpeciesCMeans$Max[7]+0.15)

ggplot(AHSingleSpeciesCMeans, aes(x=Parameter, y=Mean))+
	geom_point(data=ChenSubset,size=3,alpha=0.5,aes(x=Simple.Parameter,y=Activation.energy..eV.,colour="Chen et al. (2014)\ndataset"))+
	geom_point(data=SingleSpeciesM, size=3,alpha=0.5, aes(x=Simple.Parameter, y=Activation.energy, colour=Heterotroph.or.Autotroph))+
	geom_point(data=filter(AHSingleSpeciesCMeans,Heterotroph.or.Autotroph=="heterotroph"),size=2, aes(colour="heterotroph (mean)"))+
	geom_errorbar(data=filter(AHSingleSpeciesCMeans,Heterotroph.or.Autotroph=="heterotroph"),aes(ymin=Mean-ConfidenceInterval,ymax=Mean+ConfidenceInterval,width=0.1))+
	geom_text(data=filter(AHSingleSpeciesCMeans,Heterotroph.or.Autotroph=="heterotroph"),y=filter(AHSingleSpeciesCMeans,Heterotroph.or.Autotroph=="heterotroph")$Height,aes(label=N))+
	geom_text(data=filter(AHSingleSpeciesCMeans,Parameter=="resource growth rate"),y=max(ChenSubset$Activation.energy..eV.)+0.25,colour="turquoise4",aes(label=paste(NStrains,"Strains",sep = " ")))+
	geom_point(data=filter(AHSingleSpeciesCMeans,Heterotroph.or.Autotroph=="autotroph"),size=2,aes(colour="autotroph (mean)"))+
	geom_errorbar(data=filter(AHSingleSpeciesCMeans,Heterotroph.or.Autotroph=="autotroph"),colour="forestgreen",aes(ymin=Mean-ConfidenceInterval,ymax=Mean+ConfidenceInterval,width=0.1))+
	geom_text(data=filter(AHSingleSpeciesCMeans,Heterotroph.or.Autotroph=="autotroph"),y=filter(AHSingleSpeciesCMeans,Heterotroph.or.Autotroph=="autotroph")$Height,aes(label=filter(AHSingleSpeciesCMeans,Heterotroph.or.Autotroph=="autotroph")$N),colour="forestgreen")+
	theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
	xlab("Parameter")+ylab("Activation Energy")+ggtitle("MacArthur model parameters")+
	geom_text(data=filter(AHSingleSpeciesCMeans,Heterotroph.or.Autotroph=="heterotroph"),y=filter(AHSingleSpeciesCMeans,Heterotroph.or.Autotroph=="heterotroph")$Height,aes(label=N))+
	scale_colour_manual("",limits=c("autotroph","Chen et al. (2014)\ndataset","autotroph (mean)","heterotroph","heterotroph (mean)"),values = c("lawngreen","aquamarine3","forestgreen","darkgrey","black"))
ggsave("figures/MacArthur_with_Chen_data_v2.png",height=6,width = 10)

#Tilman Model
source("R-scripts/Tilman Functions.R")

Tilman<-read.csv("data_raw/tilman_v8.csv")
ChenSubset$Grouped.Parameter<-rep("growth rate", length(ChenSubset$Activation.energy..eV.))

#Subset of the data containing E estimates for individual species only (laboratory data)
SingleSpeciesT<-filter(Tilman,Single.species.or.multiple.=="single species"&Experiment.Type=="Lab")

#Count species numbers for each category
NkGrowth<-length(unique(filter(SingleSpeciesT,Grouped.Parameter=="k: growth")$Species))
NkUptake<-length(unique(filter(SingleSpeciesT,Grouped.Parameter=="k: nutrient uptake")$Species))
NTMortality<-length(unique(filter(SingleSpeciesT,Grouped.Parameter=="mortality rate")$Species))
NTGrowth<-length(unique(filter(SingleSpeciesT,Grouped.Parameter=="growth rate")$Species))
NTConsumption<-length(unique(filter(SingleSpeciesT,Grouped.Parameter=="consumption rate")$Species))

TSingleSpeciesCMeans<-data.frame("Parameter"=c("k: growth",
											   "k: nutrient uptake",
											   "mortality rate",
											   "growth rate",
											   "consumption rate"))
TSingleSpeciesCMeans$Mean<-listcmeans(SingleSpeciesT,ChenSubset)
TSingleSpeciesCMeans$StandardError<-listcses(SingleSpeciesT,ChenSubset)
TSingleSpeciesCMeans$ConfidenceInterval<-listcconfs(SingleSpeciesT,ChenSubset)
TSingleSpeciesCMeans$NumSpecies<-c(NkGrowth,NkUptake,NTMortality,NTGrowth,NTConsumption)
TSingleSpeciesCMeans$N<-create.n(TSingleSpeciesCMeans)
TSingleSpeciesCMeans$Max<-listcmax(SingleSpeciesT,ChenSubset)
TSingleSpeciesCMeans$Min<-listcmin(SingleSpeciesT,ChenSubset)
TSingleSpeciesCMeans$Height<-c(TSingleSpeciesCMeans$Max[1]+0.35,
							   TSingleSpeciesCMeans$Max[2]+0.35,
							   NA,
							   TSingleSpeciesCMeans$Min[4]-0.35,
							   TSingleSpeciesCMeans$Max[5]+0.35)

ggplot(TSingleSpeciesCMeans, aes(x=Parameter, y=Mean))+
	geom_point(data=ChenSubset, size=3,alpha=0.5, aes(x=Grouped.Parameter, y=Activation.energy..eV.,colour="Chen et al. (2014) \ndataset"))+
	geom_point(data=SingleSpeciesT,colour="darkgrey",size=3,alpha=0.7, aes(x=Grouped.Parameter,y=Activation.Energy))+
	geom_point(size=2)+geom_errorbar(aes(ymin=Mean-ConfidenceInterval,ymax=Mean+ConfidenceInterval,width=0.1))+
	theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
	xlab("Parameter")+ylab("Activation Energy")+ggtitle("Tilman model parameters")+
	geom_text(y=TSingleSpeciesCMeans$Height,aes(label=N))+
	geom_text(data = filter(TSingleSpeciesCMeans,Parameter=="growth rate"),y=filter(TSingleSpeciesCMeans,Parameter=="growth rate")$Max+0.35,colour="turquoise4",aes(label=paste(NStrains,"Strains",sep = " ")))+
	scale_colour_manual("",values = "aquamarine3")+ylim(min(TSingleSpeciesCMeans$Height,na.rm = TRUE),max(TSingleSpeciesCMeans$Max,na.rm=TRUE)+0.35)

ggsave("Figures/Tilman_with_Chen_data_v2.png", height=6,width=10)

#Plots showing Chen et al. (2014) dataset standard errors
PanelA<-ggplot(ChenSubset,aes(x=Standard.error.of.activation.energy,y=Activation.energy..eV.))+
	geom_point()+xlab("Standard error")+ylab("Activation energy of growth rate")+
	theme_bw()

jitter=position_jitter(width = 0.5,height = 0,seed = 555)
PanelB<-ggplot(ChenData,aes(x=1,y=Activation.energy..eV.))+
	geom_point(position = jitter,size=2,aes(colour=Standard.error.of.activation.energy))+theme_bw()+
	theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),panel.grid.major.x = element_blank())+
	scale_color_viridis("Standard Error",na.value="darkgrey")+xlab("")+ylab("Activation energy of growth rate")

arranged<-ggarrange(PanelA,PanelB,
					labels = c("A", "B"),
					ncol = 2, nrow = 1)
annotate_figure(arranged,top="Chen et al. (2014) Data")

ggsave("figures/Chen_data_SE_panels_v1.png",height=6,width=10)

PanelA+ggtitle("Chen et al. (2014) Data")+theme(plot.title = element_text(hjust=0.5))
ggsave("figures/Chen_data_SE_scatter_v1.png", height=6, width = 6)
