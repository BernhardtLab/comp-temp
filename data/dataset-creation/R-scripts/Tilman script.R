library(dplyr)
library(ggplot2)
source("R-scripts/Tilman Functions.R")

Tilman<-read.csv("data_raw/tilman_v8.csv")

#Subset of the data containing E estimates for individual species only
#(from lab experiments with >2 temperature treatments)
SingleSpeciesT<-filter(Tilman,Single.species.or.multiple.=="single species"&Experiment.Type=="Lab"&Number.of.Temps.Measured>2)

#Count species numbers for each category
NkGrowth<-length(unique(filter(SingleSpeciesT,Grouped.Parameter=="k: growth")$Species))
NkGrowth
NkUptake<-length(unique(filter(SingleSpeciesT,Grouped.Parameter=="k: nutrient uptake")$Species))
NkUptake
NTMortality<-length(unique(filter(SingleSpeciesT,Grouped.Parameter=="mortality rate")$Species))
NTMortality
NGrowth<-length(unique(filter(SingleSpeciesT,Grouped.Parameter=="growth rate")$Species))
NGrowth
NTConsumption<-length(unique(filter(SingleSpeciesT,Grouped.Parameter=="consumption rate")$Species))
NTConsumption

#Create dataframe for means
TSingleSpeciesMeans<-data.frame("Parameter"=c("k: growth",
                                              "k: nutrient uptake",
                                              "mortality rate",
                                              "growth rate",
                                              "consumption rate"))
TSingleSpeciesMeans$Mean<-listtmeans(SingleSpeciesT)
TSingleSpeciesMeans$StandardError<-listtses(SingleSpeciesT)
TSingleSpeciesMeans$ConfidenceInterval<-listtconfs(SingleSpeciesT)
TSingleSpeciesMeans$NumSpecies<-c(NkGrowth,NkUptake,NTMortality,NGrowth,NTConsumption)
TSingleSpeciesMeans$N<-create.n(TSingleSpeciesMeans)
TSingleSpeciesMeans$Max<-listtmax(SingleSpeciesT)

#Plot the data
ggplot(TSingleSpeciesMeans, aes(x=Parameter, y=Mean))+
  geom_point(data=SingleSpeciesT, size=3, alpha=0.5, colour="darkgrey", aes(x=Grouped.Parameter, y=Activation.Energy))+
  geom_point(size=2)+ geom_errorbar(aes(ymin=Mean-StandardError,ymax=Mean+StandardError,width=0.1))+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
  xlab("Parameter")+ylab("Activation Energy")+ggtitle("Tilman model parameters")+
  geom_text(y=TSingleSpeciesMeans$Max+0.15,aes(label=N))+
  ylim(min(SingleSpeciesT$Activation.Energy-0.05),max(SingleSpeciesT$Activation.Energy)+0.1)
#save the plot
ggsave("figures/Tilman_Basic_v2.png", height=6, width=10)

ggplot(TSingleSpeciesMeans, aes(x=Parameter, y=Mean))+
  geom_point(data=SingleSpeciesT, size=3, alpha=0.5, colour="darkgrey", aes(x=Grouped.Parameter, y=Activation.Energy))+
  geom_point(size=2)+ geom_errorbar(aes(ymin=Mean-ConfidenceInterval,ymax=Mean+ConfidenceInterval,width=0.1))+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
  xlab("Parameter")+ylab("Activation Energy")+ggtitle("Tilman model parameters")+
  geom_text(y=TSingleSpeciesMeans$Max+0.15,aes(label=N))+
  ylim(min(SingleSpeciesT$Activation.Energy-0.05),max(SingleSpeciesT$Activation.Energy)+0.1)
ggsave("figures/Tilman_Confidence_v2.png", width = 10,height = 6)

#Colour-code different nutrients

#Split nutrient specific data from the rest
NutrientSpecific<-filter(SingleSpeciesT,Grouped.Parameter=="consumption rate"|Grouped.Parameter=="k: nutrient uptake")
NonSpecific<-filter(SingleSpeciesT,Grouped.Parameter!="consumption rate"&Grouped.Parameter!="k: nutrient uptake")

#Plot the data
ggplot(TSingleSpeciesMeans, aes(x=Parameter, y=Mean))+
  geom_point(data=NonSpecific, size=3,alpha=0.5, colour="darkgrey", aes(x=Grouped.Parameter, y=Activation.Energy))+
  geom_point(data=NutrientSpecific, size=3,alpha=0.5, aes(x=Grouped.Parameter,y=Activation.Energy,colour=Nutrient))+
  geom_point(size=2)+ geom_errorbar(aes(ymin=Mean-StandardError,ymax=Mean+StandardError,width=0.1))+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
  xlab("Parameter")+ylab("Activation Energy")+ggtitle("Tilman model parameters")+
  geom_text(y=TSingleSpeciesMeans$Max+0.15,aes(label=N))+
  ylim(min(SingleSpeciesT$Activation.Energy-0.05),max(SingleSpeciesT$Activation.Energy)+0.1)
#save the plot
ggsave("figures/Tilman_split_nutrients.png", height=6, width=10)

#Plot of all data in original dataset (Messy)
TilmanMeans<-data.frame("Parameter"=c("k: growth",
                                      "k: nutrient uptake",
                                      "mortality rate",
                                      "growth rate",
                                      "consumption rate"))
TilmanMeans$Mean<-listtmeans(Tilman)
TilmanMeans$StandardError<-listtses(Tilman)

ggplot(Tilman, aes(x=Grouped.Parameter, y=Activation.Energy))+
  geom_point(aes(colour=Single.species.or.multiple.))+
  geom_point(data=TilmanMeans, colour="black", aes(x=Parameter, y=Mean))+
  geom_errorbar(data=TilmanMeans,aes(x=Parameter, y=Mean,ymin=Mean-StandardError,ymax=Mean+StandardError,width=0.1))+
  theme_bw()