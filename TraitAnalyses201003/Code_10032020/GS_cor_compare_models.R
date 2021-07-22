setwd("/Users/maohuang/Desktop/Kelp/Simulation_Study/SugarKelpBreeding/TraitAnalyses201003/data")

GS_cor_model<-read.csv("GS_cor_DifferentModelComparison_Plot_Level_traits_07172021.csv",sep=",",header=TRUE)

head(GS_cor_model)
#data wide to long

#Model and Dataset to not split apart
GS_cor<-reshape2::melt(GS_cor_model,id.vars=(c("Model","Dataset")))
	head(GS_cor)

# Boxplot by group

p1<-ggplot2::ggplot(data=GS_cor,aes(x=Model,y=value,fill=Dataset))+
	geom_boxplot()
boxplot(value~Model,data=GS_cor,ylab="Cross validation accuray")	
res.aov<-aov(value~Model+Dataset+Model*Dataset,data=GS_cor)
summary(res.aov)