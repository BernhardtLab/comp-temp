library(dplyr)
library(ggplot2)
source("R-scripts/MacArthur Functions.R")

MacArthur<-read.csv("data_raw/macarthur_v17.csv")

#Subset of the data containing E estimates for individual species only
#(from lab experiments with >2 temperature treatments)
SingleSpeciesM<-filter(MacArthur,Single.species.or.multiple=="single species"&Experiment.Type=="Lab"&Number.of.Temps.Measured>2)

#Count species numbers for each category
NResourceGrowth<-length(unique(filter(SingleSpeciesM,Simple.Parameter=="resource growth rate")$Species))
NResourceGrowth
NResourceK<-length(unique(filter(SingleSpeciesM,Simple.Parameter=="resource carrying capacity")$Species))
NResourceK
NMortality<-length(unique(filter(SingleSpeciesM,Simple.Parameter=="mortality rate")$Species))
NMortality
NConsumption<-length(unique(filter(SingleSpeciesM,Simple.Parameter=="consumption rate")$Species))
NConsumption
NEfficiency<-length(unique(filter(SingleSpeciesM,Simple.Parameter=="conversion efficiency")$Species))
NEfficiency

#Create dataframe for means
MSingleSpeciesMeans<-data.frame("Parameter"=c("resource growth rate",
                                              "resource carrying capacity",
                                              "mortality rate",
                                              "conversion efficiency",
                                              "consumption rate"))
MSingleSpeciesMeans$Mean<-listmmeans(SingleSpeciesM)
MSingleSpeciesMeans$StandardError<-listmses(SingleSpeciesM)
MSingleSpeciesMeans$ConfidenceInterval<-listmconfs(SingleSpeciesM)
MSingleSpeciesMeans$NumSpecies<-c(NResourceGrowth,NResourceK,NMortality,NEfficiency,NConsumption)
MSingleSpeciesMeans$N<-create.n(MSingleSpeciesMeans)
MSingleSpeciesMeans$Max<-listmmax(SingleSpeciesM)

#Plot the data
ggplot(MSingleSpeciesMeans, aes(x=Parameter, y=Mean))+
  geom_point(data=SingleSpeciesM, size=3,alpha=0.5,colour="darkgrey", aes(x=Simple.Parameter, y=Activation.energy))+
  geom_point(size=2)+ geom_errorbar(aes(ymin=Mean-StandardError,ymax=Mean+StandardError,width=0.1))+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
  xlab("Parameter")+ylab("Activation Energy")+ggtitle("MacArthur model parameters")+
  geom_text(y=MSingleSpeciesMeans$Max+0.15,aes(label=N))+
  ylim(min(SingleSpeciesM$Activation.energy-0.05),max(SingleSpeciesM$Activation.energy)+0.1)
#save the plot
ggsave("figures/MacArthur_Basic_Lab_v1.png", height=6, width = 10)

#Separate heterotroph and autotroph means

#Separate counts for autotrophs and heterotroph resources
NHeterotrophGrowth<-length(unique(filter(SingleSpeciesM,Simple.Parameter=="resource growth rate"&Heterotroph.or.Autotroph=="heterotroph")$Species))
NHeterotrophGrowth
NAutotrophGrowth<-length(unique(filter(SingleSpeciesM,Simple.Parameter=="resource growth rate"&Heterotroph.or.Autotroph=="autotroph")$Species))
NAutotrophGrowth
NHeterotrophK<-length(unique(filter(SingleSpeciesM,Simple.Parameter=="resource carrying capacity"&Heterotroph.or.Autotroph=="heterotroph")$Species))
NHeterotrophK
NAutotrophK<-length(unique(filter(SingleSpeciesM,Simple.Parameter=="resource carrying capacity"&Heterotroph.or.Autotroph=="autotroph")$Species))
NAutotrophK

#Create dataframe for means, as well as parameters for plot aesthetics (max, min, height)
AHSingleSpeciesMeans<-data.frame("Parameter"=c("resource growth rate",
                                               "resource growth rate",
                                               "resource carrying capacity",
                                               "resource carrying capacity",
                                               "mortality rate",
                                               "conversion efficiency",
                                               "consumption rate"))
AHSingleSpeciesMeans$Heterotroph.or.Autotroph<-c("heterotroph",
                                                 "autotroph",
                                                 "heterotroph",
                                                 "autotroph",
                                                 "heterotroph",
                                                 "heterotroph",
                                                 "heterotroph")
AHSingleSpeciesMeans$Mean<-listahmeans(SingleSpeciesM)
AHSingleSpeciesMeans$StandardError<-listahses(SingleSpeciesM)
AHSingleSpeciesMeans$ConfidenceInterval<-listahconfs(SingleSpeciesM)
AHSingleSpeciesMeans$NumSpecies<-c(NHeterotrophGrowth,NAutotrophGrowth,NHeterotrophK,NAutotrophK,NMortality,NEfficiency,NConsumption)
AHSingleSpeciesMeans$N<-create.n(AHSingleSpeciesMeans)
AHSingleSpeciesMeans$Max<-listahmax(SingleSpeciesM)
AHSingleSpeciesMeans$Min<-listahmin(SingleSpeciesM)
AHSingleSpeciesMeans$Height<-c(AHSingleSpeciesMeans$Max[1]+0.15,
                               AHSingleSpeciesMeans$Min[2]-0.15,
                               AHSingleSpeciesMeans$Min[3]-0.15,
                               AHSingleSpeciesMeans$Max[4]+0.15,
                               AHSingleSpeciesMeans$Max[5]+0.15,
                               AHSingleSpeciesMeans$Max[6]+0.15,
                               AHSingleSpeciesMeans$Max[7]+0.15)

#Plot the data with standard errors
ggplot(AHSingleSpeciesMeans, aes(x=Parameter, y=Mean))+
  geom_point(data=SingleSpeciesM, size=3,alpha=0.5, aes(x=Simple.Parameter, y=Activation.energy, colour=Heterotroph.or.Autotroph))+
  geom_point(data=filter(AHSingleSpeciesMeans,Heterotroph.or.Autotroph=="heterotroph"),size=2, aes(colour="heterotroph (mean)"))+
  geom_errorbar(data=filter(AHSingleSpeciesMeans,Heterotroph.or.Autotroph=="heterotroph"),aes(ymin=Mean-StandardError,ymax=Mean+StandardError,width=0.1))+
  geom_text(data=filter(AHSingleSpeciesMeans,Heterotroph.or.Autotroph=="heterotroph"),y=filter(AHSingleSpeciesMeans,Heterotroph.or.Autotroph=="heterotroph")$Height,aes(label=N))+
  geom_point(data=filter(AHSingleSpeciesMeans,Heterotroph.or.Autotroph=="autotroph"),size=2,aes(colour="autotroph (mean)"))+
  geom_errorbar(data=filter(AHSingleSpeciesMeans,Heterotroph.or.Autotroph=="autotroph"),colour="forestgreen",aes(ymin=Mean-StandardError,ymax=Mean+StandardError,width=0.1))+
  geom_text(data=filter(AHSingleSpeciesMeans,Heterotroph.or.Autotroph=="autotroph"),y=filter(AHSingleSpeciesMeans,Heterotroph.or.Autotroph=="autotroph")$Height,aes(label=N),colour="forestgreen")+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
  xlab("Parameter")+ylab("Activation Energy")+ggtitle("MacArthur model parameters")+
  geom_text(data=filter(AHSingleSpeciesMeans,Heterotroph.or.Autotroph=="heterotroph"),y=filter(AHSingleSpeciesMeans,Heterotroph.or.Autotroph=="heterotroph")$Height,aes(label=N))+
  ylim(min(SingleSpeciesM$Activation.energy-0.05),max(SingleSpeciesM$Activation.energy)+0.1)+
  scale_colour_manual(values = c("lawngreen","forestgreen","darkgrey","black"))
#save the plot
ggsave("figures/MacArthur_split_resources_lab_v1.png", height=6, width = 10)

#Plot data with confidence intervals
ggplot(AHSingleSpeciesMeans, aes(x=Parameter, y=Mean))+
  geom_point(data=SingleSpeciesM, size=3,alpha=0.5, aes(x=Simple.Parameter, y=Activation.energy, colour=Heterotroph.or.Autotroph))+
  geom_point(data=filter(AHSingleSpeciesMeans,Heterotroph.or.Autotroph=="heterotroph"),size=2, aes(colour="heterotroph (mean)"))+
  geom_errorbar(data=filter(AHSingleSpeciesMeans,Heterotroph.or.Autotroph=="heterotroph"),aes(ymin=Mean-ConfidenceInterval,ymax=Mean+ConfidenceInterval,width=0.1))+
  geom_text(data=filter(AHSingleSpeciesMeans,Heterotroph.or.Autotroph=="heterotroph"),y=filter(AHSingleSpeciesMeans,Heterotroph.or.Autotroph=="heterotroph")$Height,aes(label=N))+
  geom_point(data=filter(AHSingleSpeciesMeans,Heterotroph.or.Autotroph=="autotroph"),size=2,aes(colour="autotroph (mean)"))+
  geom_errorbar(data=filter(AHSingleSpeciesMeans,Heterotroph.or.Autotroph=="autotroph"),colour="forestgreen",aes(ymin=Mean-ConfidenceInterval,ymax=Mean+ConfidenceInterval,width=0.1))+
  geom_text(data=filter(AHSingleSpeciesMeans,Heterotroph.or.Autotroph=="autotroph"),y=filter(AHSingleSpeciesMeans,Heterotroph.or.Autotroph=="autotroph")$Height,aes(label=N),colour="forestgreen")+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
  xlab("Parameter")+ylab("Activation Energy")+ggtitle("MacArthur model parameters")+
  geom_text(data=filter(AHSingleSpeciesMeans,Heterotroph.or.Autotroph=="heterotroph"),y=filter(AHSingleSpeciesMeans,Heterotroph.or.Autotroph=="heterotroph")$Height,aes(label=N))+
  ylim(min(SingleSpeciesM$Activation.energy-0.05),max(SingleSpeciesM$Activation.energy)+0.1)+
  scale_colour_manual("",values = c("lawngreen","forestgreen","darkgrey","black"))
#save the plot
ggsave("figures/MacArthur_split_resources_confidence_v1.png",height=6,width=10)

#Visualizing how E differes between resource species consumed by the same consumer

#create data frame to store counts
ConsumerCounts<-data.frame("Species"=unique(filter(SingleSpeciesM,Simple.Parameter=="consumption rate")$Species))
ConsumerCounts$NumObservations<-list.o.counts(ConsumerCounts)
ConsumerCounts$NumResources<-list.r.counts(ConsumerCounts)

#Subset original data that contains multiple resource species for each consumer
MultipleResources<-filter(SingleSpeciesM,Species%in%filter(ConsumerCounts,NumResources>1)$Species&Simple.Parameter=="consumption rate")
MultipleResources<-MultipleResources%>%arrange(Resource.Species)
MultipleResources$ResourceID<-group_indices(MultipleResources,Resource.Species)
MultipleResources$Check<-MultipleResources$Resource.Species

#Plot the data
ggplot(MultipleResources, aes(x=Species,y=Activation.energy,group=Species))+
  geom_line(size=1.5,aes(colour=Species))+geom_point()+geom_text(nudge_x = 0.25,aes(label=ResourceID))+
  theme_bw()+theme(axis.text.x=element_blank(),axis.ticks.x = element_blank())+
  xlab("Consumer Species")+ylab("Activation Energy of Consumption Rate")+labs(colour="Consumer Species")
#save the plot
ggsave("figures/Multiple_resource_species_v2.png",height = 6,width = 10)

#With standard errors (messy)
dodge=position_dodge(0.3)
ggplot(MultipleResources, aes(x=Species,y=Activation.energy))+
  geom_line(size=1.5,aes(colour=Species,group=Species))+geom_point(aes(group=Resource.Species),position = dodge)+geom_text(nudge_x = 0.25,aes(label=ResourceID))+
  geom_errorbar(position = dodge,aes(ymin=Activation.energy-Error,ymax=Activation.energy+Error,width=0.2,group=ResourceID))+
  theme_bw()+theme(axis.text.x=element_blank(),axis.ticks.x = element_blank())+
  xlab("Consumer Species")+ylab("Activation Energy of Consumption Rate")+labs(colour="Consumer Species")

#Plot of all data in original dataset (Messy)
MacArthurMeans<-data.frame("Parameter"=c("resource growth rate",
                                         "resource carrying capacity",
                                         "mortality rate",
                                         "conversion efficiency",
                                         "consumption rate"))
MacArthurMeans$Mean<-listmmeans(MacArthur)
MacArthurMeans$StandardError<-listmses(MacArthur)

ggplot(MacArthur, aes(x=Simple.Parameter, y=Activation.energy))+
  geom_point(aes(colour=Single.species.or.multiple,shape=Experiment.Type))+
  geom_point(data=MacArthurMeans, colour="black", aes(x=Parameter, y=Mean))+
  geom_errorbar(data=MacArthurMeans,aes(x=Parameter, y=Mean,ymin=Mean-StandardError,ymax=Mean+StandardError,width=0.1))+
  theme_bw()