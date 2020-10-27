### Both Years
###################################
# Plot information analysis
###################################

library(rrBLUP)
rm(list=ls())
load("dataNHpi_withChk_3_sets_PhotoScore23.rdata")   ## Plot
load("dataNHim_withChk_3_sets_PhotoScore0123.rdata")  ## Indi

### !!!!Both
dataNHpi<-dataNHpiBoth_C  
dataNHim<-dataNHimboth_C
yr<-"Both"
  ls()
  head(dataNHpiBoth_C)

  
# ### !!!!2019
# dataNHpi<-dataNHpi19_C
# dataNHim<-dataNHim19_C
# yr<-"2019"
# 
# ### !!!!2020
# dataNHpi<-dataNHpi20_C
# dataNHim<-dataNHim20_C
# yr<-"2020"


#load(paste0("hMat_PedNH_CCmat_fndrMrkData_",yr,"_PhotoScore23.rdata"))
#load(paste0("hMat_PedNH_CCmat_fndrMrkData_Both_PhotoScore23_NoSGP.rdata"))
  load(paste0("hMat_PedNH_CCmat_fndrMrkData_Both_PhotoScore23_withSGP.rdata"))

spRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) all(vec > 0)))    ### Progeny sps, col 2 and 3 are all values
nSp <- length(spRows)
  nSp #244  #245
#This is already added
#dataNHpi$popChk <- ifelse(substr(dataNHpi$plotNo, 1, 1) == "Z", substr(dataNHpi$plotNo, 1, 2), "ES")  # Checks VS ES
  
  
  
dataNHpi$popChk <- ifelse(substr(dataNHpi$plotNo, 1, 1) == "Z", substr(dataNHpi$plotNo, 1, 2), "ES")  # Checks VS ES
  
#MHuang edit this:
# withinLoc: variable is 1 if the cross was between gametophytes sampled from
# the same location, zero otherwise
dataNHpi$withinLoc <- ifelse(as.vector(dataNHpi$femaParLoc) == as.vector(dataNHpi$maleParLoc), 1, 0) # WithinLoc is 1
dataNHpi$development[dataNHpi$development=="#N/A"] <-NA


# Experimental design variables should be factors
for (col in c( "Year", "plotNo","GrowDays","femaPar", "femaParLoc", "malePar", "maleParLoc", "PlantingDens", "block", "line","popChk")) 
  dataNHpi[,col] <- factor(dataNHpi[,col])
nlevels(dataNHpi$plotNo)

# Make data cols numeric. Enforce data is numeric
for (col in c("wetWgtPlot", "lengthPlot", "wetWgtPerM","percDryWgt",  "dryWgtPerM","densityBlades","AshFreedryWgtPerM","AshFDwPM")) 
  dataNHpi[,col] <- as.numeric(dataNHpi[,col])

# If percentDryWeigth is NA, then the plot WetWeightPerM and DryWeightperM should be set at NA, WetWeigthPerM may be just rope
dataNHpi[is.na(dataNHpi$percDryWgt), c("wetWgtPerM", "dryWgtPerM")] <- NA
keepRows <- !is.na(dataNHpi$percDryWgt)
  keepRows

# isSelf: variable is 1 if the cross was between gametophytes from the same 
# sporophyte, zero otherwise
## MHuang edit
fndrF1<-strsplit(as.character(dataNHpi$femaPar), split="-", fixed=T)
fndrM1<-strsplit(as.character(dataNHpi$malePar),split="-",fixed=T)
fndrF<-sapply(fndrF1, function(vec) paste(vec[1:3], collapse="-"))
fndrM<-sapply(fndrM1, function(vec) paste(vec[1:3],collapse="-"))

isSelf<-ifelse(fndrF==fndrM,1,0)
dataNHpi$isSelf <- isSelf 
  str(dataNHpi)

# Experimental design variables should be factors
for (col in c("Year","plotNo","GrowDays", "femaPar", "femaParLoc", "malePar", "maleParLoc", "line", "block", "date", "popChk", "withinLoc", "isSelf")) 
  dataNHpi[,col] <- as.factor(dataNHpi[,col])

### MHuang edit this
#dataNHpiES<-dataNHpi[dataNHpi$popChk=="ES",]
#dataNHpiChk<-dataNHpi[dataNHpi$popChk!="ES",]

#dataNHpi<-rbind(dataNHpiES,dataNHpiChk)  # 569 rows ES + 47 rows of Chk
### 1. Naive heritability using amatrix, Analyze GOM, BOTH years
# msX <- model.matrix( ~ line%in%Year+block%in%Year+Year+popChk, data=dataNHpi)  ### !!!! nested effects???? DID NOT give the full Ranks of X
# msX <- msX[, apply(msX, 2, function(v) !all(v == 0))]
# Construction of Z assumes that all checks are at the bottom of the dataset
# Z has dimensions nrow(dataNHpi) x ncol(aMat)
  
# the first part all 0s, diag() identity matrix
# modify msZ to add the checks
# two same cross for a different plot would be having 
### aMat has fndr+GP+ all the plot SPs in it.
### sort 
  
dataNHpi_RMchk<-dataNHpi[!dataNHpi$crossID=="Check",] # Subset without chk levels
dataNHpi_RMchk$Crosses<-as.factor(as.character(dataNHpi_RMchk$Crosses)) # This RMed the check levels
  dim(dataNHpi_RMchk)
  str(dataNHpi_RMchk)
nrowplot<-nrow(dataNHpi_RMchk)
nrowchk<-nrow(dataNHpi) - nrowplot

# ## 1. fndrs and GPs
# fndr_GPMat<-matrix(0, nrowplot, nrow(aMat) - nSp)
#   dim(fndr_GPMat)
# ## 2. plot level Sps
# msZ<-model.matrix(~ -1 +phenoNamesFact, data=dataNHpi_RMchk)
# nSpMat<-model.matrix( ~ -1+Crosses,data=dataNHpi_RMchk)  # No intercept; sorted inside the 
#   dim(nSpMat)
#   nSpMat[1:5,1:6]
# library(stringr)  
# colnames(nSpMat)<-str_replace(colnames(nSpMat),"Crosses","")  
# ### the Crosses in nSpMat,make them the same order as those crosses in aMat
# nSpMat<-nSpMat[,colnames(hMat[,(ncol(fndr_GPMat)+1):ncol(hMat)])]
# ## 4. make msZ
# msZ <- rbind(cbind(fndr_GPMat,nSpMat),chkSpMat )
#   dim(msZ)
#   dim(aMat)
#   which(msZ[19,]==1) #325
#   which(msZ[129,]==1) #325
#   hMat[325:326,325:326]
### Automate the msZ
#load("outCovComb_files.Rdata")
#identical(rownames(outCovComb1),rownames(outCovComb4)) #### 

#load("outCovComb4_09282020.Rdata")

load("outCovComb4_10012020.Rdata")

outCovComb<-outCovComb4   ###!!!!!!
  dim(outCovComb)
phenoNamesFact<-factor(dataNHpi_RMchk$Crosses,levels=rownames(outCovComb))   #colnames are sorted alphabetically for the outCovComb 
msZ0<-model.matrix(~-1 +phenoNamesFact,data=dataNHpi_RMchk)  

#### Does this matter??
identical(colnames(outCovComb4)[385:628],levels(dataNHpi_RMchk$Crosses)) ####this= FALSE !!!Ordering of crosses in levels dataNHPi differs from outCovComb1


## 3. chks rows
chkSpMat<-matrix(0, nrowchk, nrow(outCovComb))
  dim(chkSpMat)
#  
msZ<-rbind(msZ0,chkSpMat)
  dim(msZ)

# Calculate heritability from the mixed.solve output
heritability <- function(msOut){
  return(msOut$Vu / (msOut$Vu + msOut$Ve))
}

write.csv(dataNHpi,paste0("dataNHpi_Last_Used_in_Model_",yr,".csv"))

########### 6. FINALE !! RUN Model WITHOUT blade density as covariate, but use hMat as the relationship matrx

hMat<-outCovComb4 ###!!!!!!!!!aMatInitial,weight!!!!!!!!!!!!!

#apply(msX, 2, function(v) all(v == 0))
#   # ###!!!!!!!!!!!!!!!!!! Within Year, 2019, 2020 !!!!!!!!!!!!!
# 
# msX <- model.matrix( ~ line+block+popChk, data=dataNHpi)
# msX <- msX[, apply(msX, 2, function(v) !all(v == 0))]

  # AshFDwPMImpute<-ifelse(is.na(dataNHpi$AshFDwPM),0.5948*dataNHpi$dryWgtPerM-0.0139,dataNHpi$AshFDwPM)
  # dataNHpi$AshFDwPM=AshFDwPMImpute
  
  dataNHpi<-droplevels(dataNHpi) 
  
  msX1 <- model.matrix( ~ line%in%Year+block%in%Year+Year, data=dataNHpi)  ### !!!! No popChk because no checks for ashes
  msX1 <- msX1[, apply(msX1, 2, function(v) !all(v == 0))]
  
  library(Matrix)
  rankMatrix(msX1)
  qr(msX1)$rank
  
  msOutAshFDwPM<-mixed.solve(y=dataNHpi$AshFDwPM,Z=msZ, K=hMat, X=msX1, SE=T)
  msOutAshOnly<-mixed.solve(y=dataNHpi$Ash,Z=msZ,K=hMat,X=msX1,SE=T)
  
  ### This two differs between years. put in the beginning to use them
  ###!!!!!!!!!!!!!!!!!! Both Years !!!!!!!!!!!
  msX <- model.matrix( ~ line%in%Year+block%in%Year+Year+popChk, data=dataNHpi)  ### !!!! has popChk for the traits
  msX <- msX[, apply(msX, 2, function(v) !all(v == 0))]
  
  library(Matrix)
  rankMatrix(msX)
  qr(msX)$rank
  
  msOutDBh <- mixed.solve(y=dataNHpi$densityBlades, Z=msZ, K=hMat, X=msX, SE=T)
  msOutWWPh <- mixed.solve(y=log(dataNHpi$wetWgtPlot+1), Z=msZ, K=hMat, X=msX, SE=T)
  msOutDWPMh <- mixed.solve(y=log(dataNHpi$dryWgtPerM+1), Z=msZ, K=hMat, X=msX, SE=T)
  msOutPDWh <- mixed.solve(y=dataNHpi$percDryWgt, Z=msZ, K=hMat, X=msX, SE=T)
  
  h2hMat <- c(heritability(msOutAshFDwPM),heritability(msOutWWPh), heritability(msOutDWPMh), heritability(msOutPDWh),heritability(msOutDBh))
  names(h2hMat) <- c("AshFDwPM","wetWgtPlot", "dryWgtPerM", "percDryWgt","densityBladesPerM")
      round(h2hMat,3) 
      
  
allBLUPs <- cbind(msOutAshFDwPM$u,msOutAshOnly$u,msOutDWPMh$u,msOutWWPh$u,  msOutPDWh$u)
colnames(allBLUPs) <- c("AshFDwPM","AshOnly","DWpM","WWP",  "PDW")
  write.csv(allBLUPs,"allBLUPs_PlotsOnly_noSGP.csv")  
   
### Combining all the GPs that have the GEBV
SProg_GPs<-allBLUPs[grep("UCONN-S",rownames(allBLUPs)),]  ### Grep only the SProg_GPs
allBLUPs0<-allBLUPs
allBLUPs<-allBLUPs[!rownames(allBLUPs)%in%rownames(SProg_GPs),]
####

sampledSP <- which(nchar(rownames(allBLUPs)) < 11)  # The fndr SP is max of 10 characters
releaseGP <- nchar(rownames(allBLUPs))   
releaseGP <- which(10 < releaseGP & releaseGP < 15)  #FG and MG from fndr SP is max of 14 characters
progenySP <- which(nchar(rownames(allBLUPs)) > 15)   #the Cross name is larger than 15 characters


fndNames <- rownames(allBLUPs)[sampledSP]
hasDesc <- sapply(fndNames, function(n) grep(n, rownames(allBLUPs)[progenySP]))
  
nDesc <- sapply(hasDesc, length)
  nDesc
  
GP_allBLUPs<-as.data.frame(rbind(allBLUPs[releaseGP,],SProg_GPs)) ## !!! Added the SProg_GPs  
#GP_allBLUPs<-as.data.frame(rbind(allBLUPs[releaseGP,]))
SP_allBLUPs<-as.data.frame(allBLUPs[progenySP,])
SP_allBLUPs<-SP_allBLUPs[order(-SP_allBLUPs$AshFDwPM),]

#write.csv(SP_allBLUPs,paste0("GOM_RepeatCrosses_",yr,"Yr.csv")) 
#RepeatCrosses<-read.csv("GOM_RepeatCrosses_BothYr.csv",sep=",",header=TRUE)
RepeatCrosses<-SP_allBLUPs

#### Adding sugar content
# This was what's saved
#save(dataRegion19_C,dataRegion20_C,dataRegionBoth_C,file="dataRegion_withChk_3_sets_PhotoScore23.rdata")

load("dataRegion_withChk_3_sets_PhotoScore23.rdata")    ### !!!!
dataRegion<-dataRegionBoth_C
Sugars<-dataRegion[,colnames(dataRegion)%in%c("Crosses","Total.Sugars","Mannuronic","Guluronic","Glucuronic","Mannitol","Fucose" ,"Glucan","Mannan","Galactan")]
Sugars<-Sugars[, apply(Sugars, 2, function(v) !all(is.na(v)))]
  dim(Sugars)
RepeatCrosses2<-merge(RepeatCrosses,Sugars,by.x="row.names",by.y="Crosses",all.x=TRUE)
  head(RepeatCrosses2)
colnames(RepeatCrosses2)[1]<-"Crosses"  

library(expss)
GPnames<-str_split_fixed(string=as.character(RepeatCrosses2$Crosses), "x", 2)
  head(GPnames)
RepeatCrosses2$FG<-GPnames[,1]
RepeatCrosses2$MG<-GPnames[,2]
  head(RepeatCrosses2)
RepeatCrosses2<-RepeatCrosses2[order(-RepeatCrosses2$AshFDwPM),]

write.csv(RepeatCrosses2,"GOM_RepeatCrosses_AddedSugar_10012020_noSGPs.csv")



### Generate new crosses for a trait
  traitname<-c("AshFDwPM","DWpM","AshOnly")  #################!!!!!!!!!
  ListGPs<-rownames(GP_allBLUPs) 
  getGPsex <- function(gpName){
    sexNum <- strsplit(gpName, "-")[[1]][4]    ## split text and get the 4th element which is MG/FGs
    sex <- substring(sexNum, 1, 1)
    return(ifelse(sex %in% c("F", "M"), sex, NA))
  }
  gpSex <- sapply(ListGPs, getGPsex) ### apply to each element of the vector
  
  gpSex2<-matrix(nrow=length(gpSex),ncol=2)
  gpSex2[,1]<-names(gpSex)
  gpSex2[,2]<-gpSex
  gpSex2<-as.data.frame(gpSex2)
  
  colnames(gpSex2)<-c("GPList","Sex")
  gpSex2<-gpSex2[order(gpSex2$Sex),]
  
  gpF<-gpSex2[gpSex2$Sex=="F",]  #data frame
  gpM<-gpSex2[gpSex2$Sex=="M",]  #data frame
  
  CrossF<-NULL
  for (i in 1:nrow(gpF)){
    CrossFemal<-rep(as.character(gpF[i,1]),times=nrow(gpM))
    CrossF<-c(CrossF,CrossFemal)
  }
  CrossM<-rep(as.character(gpM$GPList),times=nrow(gpF))
  Crosses<-as.data.frame(cbind(CrossF,CrossM))

  for (col in c("CrossF","CrossM")){
    Crosses[,col]<-as.character(Crosses[,col])
  }
  head(Crosses)
  head(allBLUPs)
for (i in 1:length(traitname)){
  Crosses[,3*i]<-vlookup(Crosses$CrossF,dict=allBLUPs0,result_column=i,lookup_column = "row.names") # Female BLUPs
  colnames(Crosses)[3*i]<-paste0("F_",colnames(allBLUPs0)[i])
  
  Crosses[,(3*i+1)]<-vlookup(Crosses$CrossM,dict=allBLUPs0,result_column=i,lookup_column = "row.names") # Male BLUPs
  colnames(Crosses)[3*i+1]<-paste0("M_",colnames(allBLUPs0)[i])
  
  Crosses[,(3*i+2)]<-rowMeans(cbind(Crosses[,3*i],Crosses[,(3*i+1)]))
  colnames(Crosses)[3*i+2]<-paste0("Predicted_",colnames(allBLUPs0)[i])
}  
  
# Crosses$BLUPsF<-allBLUPs[Crosses$CrossF,colnames(allBLUPs)%in%traitname]
# Crosses$BLUPsM<-allBLUPs[Crosses$CrossM,colnames(allBLUPs)%in%traitname]
# Crosses$BLUPsnewSP<-rowMeans(Crosses[,colnames(Crosses)%in%c("BLUPsF","BLUPsM")])

Crosses$FGLoc<-str_split_fixed(string=Crosses$CrossF,"-",4)[,2]
Crosses$MGLoc<-str_split_fixed(string=Crosses$CrossM,"-",4)[,2]

SNELoc<-c("TI", "FW", "PI", "FI", "BL", "DB", "LR")
SNERegion<-cbind(rep("SNE",times=length(SNELoc)),SNELoc)
colnames(SNERegion)<-c("Region","Loc")
  head(SNERegion)
  str(SNERegion)
Crosses$FRegion<-vlookup(Crosses$FGLoc,dict=SNERegion,result_column ="Region",lookup_column = "Loc")
Crosses$MRegion<-vlookup(Crosses$MGLoc,dict=SNERegion,result_column ="Region",lookup_column = "Loc")
  head(Crosses)
Crosses<-Crosses[order(-Crosses$Predicted_DWpM),]
  head(Crosses)
write.csv(Crosses,paste0("GOM_NewCrosses_AddedTraits_",yr,"Yr_10012020_noSGPs.csv"))     #### !!!!!!!!!!!!!
  
  

##### Update checking on biomass
#### Checking on biomass
summary<-read.csv("Active_culture_biomass_assesment_for_Mao_100120.csv",sep=",",header=T,na.strings=c("","NA"))  ###!!!!
  dim(summary)  
  head(summary)
summary$GametophyteID<-as.character(summary$GametophyteID)
  str(summary)
  
RepeatCrosses<-read.csv("GOM_RepeatCrosses_AddedSugar_10012020.csv",sep=",",header=T,row.names=1)

RepeatCrosses2$FGBio<-vlookup(RepeatCrosses2$FG,dict=summary,result_column="X..Crosses",lookup_column ="GametophyteID")
RepeatCrosses2$MGBio<-vlookup(RepeatCrosses2$MG,dict=summary,result_column="X..Crosses",lookup_column ="GametophyteID")

RepeatCrosses2$Manual<-ifelse(is.na(RepeatCrosses2$FGBio) |is.na(RepeatCrosses2$MGBio) | RepeatCrosses2$FGBio==0 | RepeatCrosses2$MGBio==0,"RM","Keep")


Crosses<-read.csv("GOM_NewCrosses_AddedTraits_BothYr_10012020_AddBiomass.csv",sep=",",header=T,row.names=1) 
Crosses$FBio<-vlookup(Crosses$CrossF,dict=summary,result_column = "X..Crosses",lookup_column="GametophyteID")
Crosses$MBio<-vlookup(Crosses$CrossM,dict=summary,result_column = "X..Crosses",lookup_column="GametophyteID")
Crosses$Manual<-ifelse(is.na(Crosses$FBio) |is.na(Crosses$MBio) | Crosses$FBio==0 | Crosses$MBio==0,"RM","Keep")


write.csv(RepeatCrosses2,"GOM_RepeatCrosses_AddedSugar_10012020_AddBiomass.csv")
write.csv(Crosses,paste0("GOM_NewCrosses_AddedTraits_BothYr_10012020_AddBiomass.csv"))     #### !!!!!!!!!!!!!

  head(Crosses)

  


############# Individual
for (col in c( "Year", "plotNo")) 
  dataNHim[,col] <- factor(dataNHim[,col])

tail(dataNHim)
dim(dataNHim)
tail(dataNHpi)

for (col in c("bladeLength", "bladeMaxWidth", "bladeThickness", "stipeLength", "stipeDiameter"))  
  dataNHim[,col] <- as.numeric(dataNHim[,col])

#### dataNHim of plotNo needs to be the same as that in dataNHpi, They were sorted alphabetically before saving to the rdata object
################## the order of plot in msXim/dataNHim matter $$$$

##Complicated by the fact that there are no data for some im plots
#emZpl <- model.matrix( ~ -1 + plotNo, data=dataNHim)  #### Does this need to be connecting to Crosses???
#rownames(dataNHim)<-paste0("im_",1:nrow(dataNHim))
#rownames(emZpl)<-rownames(dataNHim)
# levels(dataNHim$plotNo) <- levels(dataNHpi$plotNo)
##OR USE THIS:  dataNHim$plotNo <- factor(dataNHim$plotNo, levels=levels(dataNHpi$plotNo))
########## Now the levels in dataNHim is more than those in dataNHpi, because the dataNHpi is only using plots with photo score 2,3. But dataNHim is having all plots that produced individual measurements

dataNHim2<-subset(dataNHim,dataNHim$plotNo%in%dataNHpi$plotNo) #!!!! Each plotNo is unique
dim(dataNHim2)
dim(dataNHim)
dataNHim2$plotNo<-factor(dataNHim2$plotNo)
nlevels(dataNHim2$plotNo) 
nlevels(dataNHpi$plotNo)
nlevels(dataNHim$plotNo)

levels(dataNHim2$plotNo) <- levels(dataNHpi$plotNo)
dataNHim0<-dataNHim
dataNHim<-dataNHim2

### No Need to be  changed to Crosses, individual connects with plots
msXim <- model.matrix( ~ -1 + plotNo, data=dataNHim)  # nrow=nrow(dataNHim), ncol=number of plotNo=nlevels(plotNo)

#rownames(msXim)<-rownames(dataNHim)
dim(msXim)
dim(dataNHim)
# row(msXim) is the individual measurements, ordered alphabeticaly in dataNHim; col(msXim) is the levels of plotNo matched to dataNHpi$plotNo
# row(msZ) is the plotsNo, ordered alphabetically in dataNHpi; col(msZ) is the order of fndr,FG,MG,plots(SP progeny) in vector u, and in hMat
# row(msXdb) is the plotsNo in dataNHpi,col(msXdb) are the different factors
# rownames(dataNHpi)<-dataNHpi$plotNo  
# msXdb <- model.matrix( ~ densityBlades + line + block + popChk, data=dataNHpi) ### adding densityBlades, msXdb
# 
#msX <- model.matrix( ~ line%in%Year+block%in%Year+Year+popChk, data=dataNHpi)
msZim <- emZsp <- msXim %*% msZ
#msXim <- msXim %*% msXdb
msXim<-msXim%*%msX

dim(msX)
dim(msXim)
dim(msZim)
dim(dataNHim)

msXim<-msXim[, apply(msXim, 2, function(v) !all(v == 0))]

msOutBLh <- mixed.solve(y=log(dataNHim$bladeLength+1), Z=msZim, K=hMat, X=msXim, SE=T)
msOutBMWh <- mixed.solve(y=log(dataNHim$bladeMaxWidth+1), Z=msZim, K=hMat, X=msXim, SE=T)
msOutBTh <- mixed.solve(y=log(dataNHim$bladeThickness+1), Z=msZim, K=hMat, X=msXim, SE=T)
msOutSLh <- mixed.solve(y=log(dataNHim$stipeLength+1), Z=msZim, K=hMat, X=msXim, SE=T)
msOutSDh <- mixed.solve(y=log(dataNHim$stipeDiameter+1), Z=msZim, K=hMat, X=msXim, SE=T)
h2hMatim <- c(heritability(msOutBLh), heritability(msOutBMWh), heritability(msOutBTh), heritability(msOutSLh), heritability(msOutSDh)) #
names(h2hMatim) <- c("bladeLength", "bladeMaxWidth",  "bladeThickness", "stipeLength", "stipeDiameter") 

round(h2hMatim,3)

names(h2hMat) 
names(h2hMatim) <- c("bladeLength", "bladeMaxWidth", "bladeThickness", "stipeLength", "stipeDiameter")

allh2h<-c(h2hMat,h2hMatim)
allh2h
write.csv(allh2h,paste0("heritability_all_",yr,".csv"))


#rownames(allBLUPs) is the same as that in u, and hMat
allBLUPs <- cbind(msOutAshFDwPM$u,msOutAshOnly$u,msOutDWPMh$u,msOutWWPh$u,  msOutPDWh$u, msOutDBh$u, msOutBLh$u, msOutBMWh$u,  msOutBTh$u, msOutSLh$u, msOutSDh$u)
colnames(allBLUPs) <- c("AshFDwPM","AshOnly","DWpM","WWP",  "PDW", "BDns", "BLen", "BMax",  "BThk", "SLen", "SDia")

pdf(paste0("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/Phenotypic_Analysis/","CorrPlot_",yr,"_PhotoScore23.pdf"))
corrplot::corrplot.mixed(cor(allBLUPs), diag="n", tl.cex=0.9, tl.col=1)
dev.off()



  # Different colors for fndr, GP, progeny SP crosses
colSPGP <- c(rep("black", length(sampledSP)), rep("dark green", length(releaseGP)), rep("dark red", length(progenySP)))
colSPGP[which(nDesc == 0)] <- "grey"
  
pdf(paste0("DWpMvsIndividual",yr,"_PhotoScore23.pdf"))
plot(allBLUPs[,2], pch=16, col=colSPGP, xlab="Individual number", ylab="Dry weight per plot") # !!The second col is DWpM
  dev.off()
  
rnPos <- rownames(allBLUPs)[intersect(which(nDesc == 0), which(allBLUPs[,2]>0))] #!! The second col is DWpM
rnNeg <- rownames(allBLUPs)[intersect(which(nDesc == 0), which(allBLUPs[,2]<0))] #!! The second col is DWpM
locPos <- table(substring(rnPos, 6, 7))
locNeg <- table(substring(rnNeg, 6, 7))
  
allBLUPsDF <- as.data.frame(allBLUPs)
allBLUPsDF <- cbind(pedigree=rownames(allBLUPs), allBLUPsDF)
  #rownames(allBLUPsDF)[progenySP] <- paste("SL",yr,"-UCONN-S", dataNHpi$plotNo[1:length(progenySP)], sep="")  #### Add the plot No to it
  ### THIS DOES NOT WORK because it has "." rathe than "-", so later it is hard to strsplit()
  # rownames(allBLUPsDF)[progenySP]<-c(dataNHpi$crossID[1:length(progenySP)])
  head(allBLUPsDF)
  tail(allBLUPsDF)
  
pdf(paste0("DWpMvsPDWsporophytes",yr,"_PhotoScore23.pdf"))
plot(allBLUPsDF$PDW[progenySP], allBLUPsDF$DWpM[progenySP], pch=16, cex=0.8, ylab="Dry weight per meter", xlab="Percent dry weight", main="Sporophytes", cex.lab=1.3)
  dev.off()
pdf(paste0("DWpMvsPDWgametophytes",yr,"_PhotoScore23.pdf"))
plot(allBLUPsDF$PDW[releaseGP], allBLUPsDF$DWpM[releaseGP], pch=16, cex=0.8, ylab="Dry weight per meter", xlab="Percent dry weight", main="Gametophytes", cex.lab=1.3)
  dev.off()
  
allBLUPsDF <- cbind(pedigree=allBLUPsDF$pedigree, index=2*allBLUPsDF$DWpM + allBLUPsDF$PDW, allBLUPsDF[,-1]) #Re
  
saveRDS(allBLUPsDF, file=paste0("allBLUPsDF_analyzeNH_",yr,"_PhotoScore23.rds"))
write.csv(allBLUPsDF,paste0("allBLUPs_DF_",yr,".csv"))
  



#########!!!!!!!!!!!!! Run this to get the outCovComb4 Only Once !!!!!!!!
source("is.square.matrix.R")
source("is.positive.definite.R")
source("is.symmetric.matrix.R")
load("GPsAmat_NA0.8_P1P2P3_09282020.Rdata")

#make positive definite: decimal points differ off-diagonal causes it being non-symmetric
# Further reduced the "SL18-LD-13-Mg-3" #GPsA<-GPsA[!rownames(GPsA)=="SL18-LD-13-Mg-3",!colnames(GPsA)=="SL18-LD-13-Mg-3"]
which(colnames(GPsA)=="SL18-PI-1-FG1")
IndiRM<-c("SL18-PI-1-FG1","3","SL18-CT1-FG3","SL18-CT1-MG2","SL18-CT1-MG3","SL18-OI-15-Female","SL18-ME-1-FG1","SL18-ME-1-MG1","SL18-ME-1-MG2")

GPsA<-GPsA[!rownames(GPsA)%in%IndiRM,!colnames(GPsA)%in%IndiRM] # 269 by 269
dim(GPsA)

fndrsA<-round(mrkRelMat)     ### Added "shrink=TRUE" when estimating the mrkRelMat2
diag(fndrsA) <- diag(fndrsA) + 1e-5
is.positive.definite(fndrsA)   ### RelationshipMatrix_1. fnders A

GPsA<-round(GPsA,digits=5)   ### RelationshipMatrix_2.GPs A
diag(GPsA) <- diag(GPsA) + 1e-5
is.positive.definite(GPsA)  

#diag(mrkRelMat2) <- diag(mrkRelMat2) + 1e-5
diag(aMat) <- diag(aMat) + 1e-5
is.positive.definite(aMat)     ### RelationshipMatrix_3: pedigree based for ALL

#save(fndrsA,GPsA,aMat,file="CovList_3_As.Rdata")
#### Run these in terminal too !!!!!
#load("CovList_3_As.Rdata")
#install.packages("CovCombR")

library(CovCombR)

#### III. Make the combined RelationshipMatrix
### 1. NO initial, 3-list
CovList<-NULL
CovList[[1]]<-fndrsA ## fndrsA
CovList[[2]]<-GPsA  ## Further RMed the "SL18-LD-13-Mg-3"
CovList[[3]]<-aMat   

### outCovComb1<-CovComb(CovList,nu=1500)   

### 2. aMat initial, 2-list ------- Not correct, Only 239 invididuals
# CovList2<-NULL
# CovList2[[1]]<-fndrsA ## fndrsA
# CovList2[[2]]<-GPsA
# outCovComb2<-CovComb(CovList2,nu=1500,Kinit=aMat) 
####

### 3. aMat initial, 3-list
### outCovComb3<-CovComb(CovList,nu=1500,Kinit=aMat) 

### 4. aMat initial, 3-list, add weight
weights<-c(2,2,1)
outCovComb4<-CovComb(CovList,nu=1500,w=weights,Kinit=aMat) 

save(outCovComb4,file="outCovComb4_10012020.Rdata")
#########!!!!!!!!!!!!! Run the above Only Once !!!!!!!!


  
##### CV !!!!!!  Done only once 
#   #This is for setting CV, with a 20 folds scheme
#   #nreps<-round(nrow(dataNHpiYr)/20)+1
#   #sets<-rep(1:20,nreps)[-1:-(nreps*20-nrow(dataNHpiYr))]
#   #sets<-sets[order(runif(nrow(dataNHpiYr)))]
## Test pop: GPs used for making crosses (repeats + new combinations) and the ones genotyped
## Adding these union sets into a matrix


  ### random sampling to make 30 NAs on the dataNHpi
  ### test<-dataNHpi with NAs
  ### test rows of the u
  ### cor(test_u and the test_dataNHpi)
  
  cycle=200
  nSample=20
  traitname<-"dryWgtPerM"
  
  hMat <- calcHmatrix(mrkRelMat, aMat, aMatFounders=rownames(mrkRelMat))
  hMat<-hMat[rownames(outCovComb4),colnames(outCovComb4)]
  identical(colnames(hMat),colnames(outCovComb4))
  aMat<-aMat[rownames(outCovComb4),colnames(outCovComb4)]
  
  ### cor=1 for outCovComb1 and outCovComb3 
  
  CovCombList<-list(aMat,hMat,outCovComb1,outCovComb3,outCovComb4)
  
  #corMeans<-matrix(nrow=length(CovCombList),ncol=2)
  
  
  cor.r<-matrix(nrow=cycle,ncol=length(CovCombList))
  colnames(cor.r)<-c("aMat","hMat","NoInit","aMatInit","aMatInit_withWeight")
  rownames(cor.r)<-c(paste0("cor.r_",1:cycle))
  
  cor.p<-matrix(nrow=cycle,ncol=length(CovCombList))
  colnames(cor.p)<-c("aMat","hMat","NoInit","aMatInit","aMatInit_withWeight")
  rownames(cor.p)<-c(paste0("cor.p_",1:cycle))
  
  #corList<-NULL
  
  
  for (i in 1:cycle){
    sampleYr<-sample(c(2019,2020),1,replace=TRUE)
    dataNHpiYr<-dataNHpi_RMchk[dataNHpi_RMchk$Year==sampleYr,]
    
    subsets<-order(runif(nrow(dataNHpiYr)))
    #print(subsets)
    pred.set<-dataNHpiYr[subsets<=nSample,]
    pred.set$Crosses<-as.factor(pred.set$Crosses)
    #print(str(pred.set[,1:15]))
    
    Y<-dataNHpi
    Y[rownames(Y)%in%rownames(pred.set),colnames(Y)==traitname]=NA
    
    for (j in 1:length(CovCombList)){
      hMat<-CovCombList[[j]]
      
      msOutDWPMh <- mixed.solve(y=log(Y[,traitname]+1), Z=msZ, K=hMat, X=msX, SE=T)
      AllhMatBLUPs<-msOutDWPMh$u
      pred.set.u<-AllhMatBLUPs[as.character(pred.set$Crosses)]
      pred.set.y<-pred.set[,colnames(pred.set)==traitname]
      
      cor.r[i,j]<-cor.test(pred.set.u,pred.set.y)$estimate
      cor.p[i,j]<-cor.test(pred.set.u,pred.set.y)$p.value
    }
  }
  
  #corList[[j]]<-cor
  corMeans<-matrix(nrow=ncol(cor.r),ncol=2)
  rownames(corMeans)<-colnames(cor.r)
  colnames(corMeans)<-c("r","p")
  
  corMeans[,1]<-t(colMeans(cor.r))
  corMeans[,2]<-t(colMeans(cor.p))
  corMeans
  
  #r   p.value (20 cycles) (same sample set)
  #aMat                0.2711160 0.3381301
  #hMat                0.3094897 0.2670548
  #NoInit              0.3030436 0.2853737
  #aMatInit            0.3030436 0.2853737
  #aMatInit_withWeight 0.3085474 0.2757711
  
  #r   p.value (200 cycles) (same sample set)
  #aMat                0.2295683 0.3912472
  #hMat                0.2576408 0.3626239
  #NoInit              0.2522064 0.3600591
  #aMatInit            0.2522064 0.3600591
  #aMatInit_withWeight 0.2544694 0.3559478
  
  #r   p.value (200 cycles) (not the same sample set)
  #NoInit              0.23319457 0.3716543
  #aMatInit            0.24153700 0.3679468
  #aMatInit_withWeight 0.23943776 0.3571202
##### CV ends here !!!!!
  
  
  
  # Best sporophytes
  nTop<-20  ### Pick the top 20
  bestSP_DWpM <- allBLUPsDF[progenySP,][order(allBLUPsDF$DWpM[progenySP], decreasing=T)[1:nTop],] # based on DWpM
  bestSP_DWpM[,-(1:2)] <- bestSP_DWpM[,-(1:2)] %>% round(3) # Only keep the numeric part
  write.table(bestSP_DWpM, paste0("BestSPbyDryWgtPerM_",yr,"_PhotoScore23.txt"), quote=F, row.names=T, col.names=T)
  
  bestSP_idx <- allBLUPsDF[progenySP,][order(allBLUPsDF$index[progenySP], decreasing=T)[1:nTop],]
  bestSP_idx[,-(1:2)] <- bestSP_idx[,-(1:2)] %>% round(3)
  write.table(bestSP_idx, paste0("BestSPbyDWpMandPDW_",yr,"_PhotoScore23.txt"), quote=F, row.names=T, col.names=T)
  
  # Best gametophytes
  bestGP_DWpM <- allBLUPsDF[releaseGP,][order(allBLUPsDF$DWpM[releaseGP], decreasing=T)[1:nTop],]
  bestGP_DWpM[,-(1:2)] <- bestGP_DWpM[,-(1:2)] %>% round(3)
  write.table(bestGP_DWpM[,-1], paste0("BestGPbyDryWgtPerM_",yr,"_PhotoScore23.txt"), quote=F, row.names=T, col.names=T)
  bestGP_idx <- allBLUPsDF[releaseGP,][order(allBLUPsDF$index[releaseGP], decreasing=T)[1:nTop],]
  bestGP_idx[,-(1:2)] <- bestGP_idx[,-(1:2)] %>% round(3)
  write.table(bestGP_idx[,-1], paste0("BestGPbyDWpMandPDW_",yr,"_PhotoScore23.txt"), quote=F, row.names=T, col.names=T)
  
  
