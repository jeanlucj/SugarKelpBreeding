### Run an ANOVA

rm(list=ls())
### The dataNHpi needs to have all the checks at the bottom in order to construct the Z matrix properly

WD<-"/Users/maohuang/Desktop/Kelp/Simulation_Study/SugarKelpBreeding/TraitAnalyses201003"

datafdr<-paste0(WD,"/data/")
setwd(WD)
  getwd()
  ls()

library(rrBLUP)

load(paste0(datafdr,"dataNHpi_withChk_3Yrs_PhotoScore123_07152021.rdata"))  
load(paste0(datafdr,"dataNHim_withChk_3Yrs_PhotoScore0123_07152021.rdata"))
  


### !!!!!!!!!! Here needs to run 3 times, each for each scenario
### !!!!Both
dataNHpi<-dataNHpi3yrs_C  
dataNHim<-dataNHim3yrs_C
yr<-"Three"
  ls()
  str(dataNHpi)
  str(dataNHim)
bothYr<-TRUE
####### These below are the data formatting procedure 
dataNHpi$densityBlades<-ifelse(dataNHpi$densityBlades==0,NA,dataNHpi$densityBlades)  # densityblades as 0 then NA
dataNHpi$popChk <- ifelse(substr(dataNHpi$plotNo, 1, 1) == "Z", substr(dataNHpi$plotNo, 1, 2), "ES")  # Checks VS ES

dataNHpi$withinLoc <- ifelse(as.vector(dataNHpi$femaParLoc) == as.vector(dataNHpi$maleParLoc), 1, 0) # WithinLoc is 1
dataNHpi$development[dataNHpi$development=="#N/A"] <-NA
  str(dataNHpi)

# Experimental design variables should be factors
for (col in c( "Year", "plotNo","Region","femaPar", "femaParLoc", "malePar", "maleParLoc", "block", "line","popChk"))  #"GrowDays"
  dataNHpi[,col] <- factor(dataNHpi[,col])
  nlevels(dataNHpi$plotNo)
  nlevels(dataNHpi$Year)
  unique(dataNHpi$Year)

# If percentDryWeigth is NA, then the plot WetWeightPerM and DryWeightperM should be set at NA, WetWeigthPerM may be just rope
dataNHpi[is.na(dataNHpi$percDryWgt), c("dryWgtPerM")] <- NA
keepRows <- !is.na(dataNHpi$percDryWgt)
  sum(keepRows)  # 385

fndrF1<-strsplit(as.character(dataNHpi$femaPar), split="-", fixed=T)
fndrM1<-strsplit(as.character(dataNHpi$malePar),split="-",fixed=T)
fndrF<-sapply(fndrF1, function(vec) paste(vec[1:3], collapse="-"))
fndrM<-sapply(fndrM1, function(vec) paste(vec[1:3],collapse="-"))

isSelf<-ifelse(fndrF==fndrM,1,0)
dataNHpi$isSelf <- isSelf 
  str(dataNHpi)

# Experimental design variables should be factors
  #"GrowDays", "date",
for (col in c("Year","plotNo","femaPar", "femaParLoc", "malePar", "maleParLoc", "line", "block",  "popChk", "withinLoc", "isSelf")) 
  dataNHpi[,col] <- as.factor(dataNHpi[,col])

####### These above are the data formatting procedure  


data<-droplevels(dataNHpi[!dataNHpi$PhotoScore==1,])
	dim(data)
	data[250:284,]
	
exptlSP <- as.character(data$popChk) == "ES"
data$entry <- as.character(data$popChk)
data$entry[exptlSP] <-as.character(data$Crosses[exptlSP])

data$group <- as.factor(ifelse(exptlSP, 1, 0))
data0<-data  


YR<-c(2021)
WithinYr<-TRUE

YR<-c(2019,2020)
WithinYr<-FALSE

YR<-c(2020,2021)
WithinYr<-FALSE

YR<-c(2019,2020,2021)
WithinYr<-FALSE


#"wetWgtPerM"    "percDryWgt"    "dryWgtPerM" "densityBlades"

data<-droplevels(data0[data0$Year%in%YR,])
	dim(data)
	str(data)
	head(data)
data$Trait<-data$wetWgtPerM ##!!!	

if(WithinYr==TRUE){
	if (yr==2021){
fitAug <- lme4::lmer(Trait ~ line*block + (1|entry:group), data=data)		
	}else {
fitAug <- lme4::lmer(Trait ~ popChk + line*block + (1|entry:group), data=data)	
		}
}else if (WithinYr==FALSE){
fitAug <- lme4::lmer(Trait ~ popChk + line*block + (1|entry:group)+(1|Year), data=data)
}

 print(aov <- anova(fitAug))
  
 lmerTest::ranova(fitAug)
 #DWpM
#2019+2020. line p-value=0.0165  entry and year both *
#2019, nothing sig
#2020, not converged
#2021, entry sig
#3years entry and year both *
#2020,2021, entry and year both *
  
  #pDW,DB,WWpM
#3years: year *  

 #### 
 WD<-"/Users/maohuang/Desktop/Kelp/Simulation_Study/SugarKelpBreeding/TraitAnalyses201003/"
#WD<-"/local/workdir/mh865/GCA_SCA/OneTime192021/"  # run in terminal

datafdr<-paste0(WD,"data/")

  pop<-"Yr19to21"

 load(paste0(datafdr,"dBLUPs_3_yrs_cal_WithinYr_",pop,".Rdata"))
##### 
Y<-WithinYr_Both_dBLUPs

Y$Year<-as.factor(WithinYr_Both_dBLUPs$Year.x)
Y$Crosses<-WithinYr_Both_dBLUPs$Crosses.x
Y$Trait<-Y$bladeThickness  ## !!!!

library(ggplot2)
plot<-ggplot(data=Y,aes(Trait,Year))+
  geom_point(aes(color=as.factor(Year)))+ 
  geom_line(aes(group=as.factor(Crosses)))
print(plot)
 