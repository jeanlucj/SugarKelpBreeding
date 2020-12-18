# GOM & SNE

rm(list=ls())
dataAll <-read.csv("dataNHpi_ManuallyAddGrowth_9.29.20_Ireland_Data_AshFDWpM.csv",sep=",",header=TRUE) ##!!!
dataAll$AshFDwPM<-(dataAll$wetWgtPerM*dataAll$percDryWgt/100)*(1-(dataAll$Ash/100))

dataAll$popChk <- ifelse(substr(dataAll$plotNo, 1, 1) == "Z", substr(dataAll$plotNo, 1, 2), "ES")  # Checks VS ES

exptlSP <- as.character(dataAll$popChk) == "ES"
dataAll$entry <- as.character(dataAll$popChk)
dataAll$entry[exptlSP] <- as.character(dataAll$plotNo[exptlSP])
dataAll$group <- as.factor(ifelse(exptlSP, 1, 0))

  ls()
  head(dataAll)
  dataAll[1:10,]
  tail(dataAll)
  dim(dataAll)
  str(dataAll)
#dataRegion<-dataAll[dataAll$Region=="SNE",] ## !!! RM SNE, 530 rows\

dataRegion<-dataAll[dataAll$Region=="GOM",] ## !!! RM SNE, 530 rows\

dataRegion <- dataRegion[order(dataRegion$plotNo),]  #### Plots in alphabetic order

dataRegion<-dataRegion[!dataRegion$crossID=="Buffer",]  ## !!! RMed Buffer lines

#dataRegion<-dataRegion[!dataRegion$PhotoScore==0,] ## !!! RM PhotoScore=0,  447 rows
dataRegion<-dataRegion[dataRegion$PhotoScore>1,] ## !!! RM PhotoScore<=1
  dim(dataRegion)
  colnames(dataRegion)
dataRegion19<-dataRegion[dataRegion$Year==2019,] ## 2019 data   !!!!!!
dataRegion19_C<-dataRegion19
dataRegion19_C<-dataRegion19_C[order(dataRegion19_C$plotNo),]  ## Order plotNo alphabetically

dataRegion20<-dataRegion[dataRegion$Year==2020,]   #
dataRegion20_C<-dataRegion20
dataRegion20_C<-dataRegion20_C[order(dataRegion20_C$plotNo),]  ## Order plotNo alphabetically

dataRegionBoth_C<-dataRegion  ## Order plotNo alphabetically
dataRegionBoth_C<-dataRegionBoth_C[order(dataRegionBoth_C$plotNo),] ## Order plotNo alphabetically

#save(dataRegion19_C,dataRegion20_C,dataRegionBoth_C,file="dataRegion_withChk_3_sets_PhotoScore23.rdata")



  dim(dataRegion19_C)
  dim(dataRegion20_C)
  dim(dataRegionBoth_C)  
  head(dataRegion20_C)

rm(list=ls())
load("dataRegion_withChk_3_sets_PhotoScore23.rdata")
  # Experimental design variables should be factors
  
  for (col in c( "Year", "plotNo","GrowDays","femaPar", "femaParLoc", "malePar", "maleParLoc", "PlantingDens", "block", "line","popChk")) 
    dataRegion[,col] <- factor(dataRegion[,col])
  nlevels(dataRegion$plotNo)
  
dataRegion$withinLoc <- ifelse(as.vector(dataRegion$femaParLoc) == as.vector(dataRegion$maleParLoc), 1, 0) # WithinLoc is 1
  dataRegion$development[dataRegion$development=="#N/A"] <-NA
  
  # Make data cols numeric. Enforce data is numeric
  for (col in c("wetWgtPlot", "lengthPlot", "wetWgtPerM","percDryWgt",  "dryWgtPerM","densityBlades","AshFDwPM")) 
    dataRegion[,col] <- as.numeric(dataRegion[,col])
  
  # If percentDryWeigth is NA, then the plot WetWeightPerM and DryWeightperM should be set at NA, WetWeigthPerM may be just rope
  dataRegion[is.na(dataRegion$percDryWgt), c("wetWgtPerM", "dryWgtPerM")] <- NA
  keepRows <- !is.na(dataRegion$percDryWgt)
  keepRows
  
  # isSelf: variable is 1 if the cross was between gametophytes from the same 
  # sporophyte, zero otherwise
  ## MHuang edit
  fndrF1<-strsplit(as.character(dataRegion$femaPar), split="-", fixed=T)
  fndrM1<-strsplit(as.character(dataRegion$malePar),split="-",fixed=T)
  fndrF<-sapply(fndrF1, function(vec) paste(vec[1:3], collapse="-"))
  fndrM<-sapply(fndrM1, function(vec) paste(vec[1:3],collapse="-"))
  
  isSelf<-ifelse(fndrF==fndrM,1,0)
  dataRegion$isSelf <- isSelf 
    str(dataRegion)
  
  # Experimental design variables should be factors
  for (col in c("Year","plotNo","GrowDays", "femaPar", "femaParLoc", "malePar", "maleParLoc", "line", "block", "date", "popChk", "withinLoc", "isSelf","entry","group")) 
    dataRegion[,col] <- as.factor(dataRegion[,col])

    write.csv(dataRegion,paste0("dataRegion_Last_Used_in_Model_SugarsBothYr.csv"))
    
   
dataRegion$Trait<-dataRegion$Total.Sugars ## !

#  DId not work
#Error: number of levels of each grouping factor must be < number of observations (problems: entry:group)
# library(lme4)
# LSMeans<-lmer(Trait~line:Year+block:Year+Year+group+(1|group/entry),data=dataRegion)
# 
# print(aov<-anova(LSMeans))


# For SNE Only
#LSMeans<-lmer(DwPM~line*block+(1|entry),data=dataRegion,control=lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"))
#Because the number of groups for checks is only 1 observation, no nested effects
# #install.packages("emmeans")
# library(emmeans)
# LSMeans<-lm(Trait~line*block+entry,data=dataRegion)
#   summary(LSMeans)
# DwPM<-lsmeans(LSMeans,"entry")  


#install.packages("emmeans")
library(emmeans)
LSMeans<-lm(Trait~line*Year+block*Year+entry,data=dataRegion)
  summary(LSMeans)
Trait<-lsmeans(LSMeans,"entry")  #entry %in% (line*Year*block)

write.csv(Trait,"dataRegionTrait_LSMeans.csv")
Trait<-read.csv("dataRegionTrait_LSMeans.csv",sep=",",header=T,row.names=1)
  head(Trait)
  Trait
library(plyr)
# Note: there are two checks can be macthed up to this in the dataRegion!!!
Trait2<-merge(Trait[,colnames(Trait)%in%c("entry","line","block","lsmean")],dataRegion[,colnames(dataRegion)%in%c("Trait","Crosses","entry")],by="entry")   

Trait2<-Trait2[order(-Trait2$lsmean),]

library(stringr)
GPName<-str_split_fixed(string=as.character(Trait2$Crosses), "x", 2)
  GPName
Trait2$FG<-GPName[,1]
Trait2$MG<-GPName[,2]
  head(Trait2)
write.csv(Trait2,"RepeatCrosses_GOM_SugarTrait.csv")


#### It turns out these lsmeans are not able to well be estimated, 
#So just use the row means


