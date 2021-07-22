###ASreml-R trait genetic correlations between years
# save(dataNHpi,All_grm,Trait_grm,file="Asreml_OtherTraits.Rdata")

###############################################  
rm(list=ls())
#WD<-"/Users/maohuang/Desktop/Kelp/GCA_SCA/" # local
#datafdr<-"/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/data/"
datafdr<-paste0("/Users/maohuang/Desktop/Kelp/Simulation_Study/SugarKelpBreeding/TraitAnalyses201003/data/")

WD<-"/local/workdir/mh865/GCA_SCA/"
datafdr<-"/local/workdir/mh865/ASreml/"

###1. Input outCovComb4_dipOrder
load(paste0(datafdr,"outCovComb4_Mix_Conden_0712_2021.Rdata"))

All_grm<-outCovComb4_MixOrder  # 950

###2. Input dataNH

#### Plot level
# load(paste0(datafdr,"dataNHpi_withChk_3_sets_PhotoScore23_UpdateAsh_0309_2021.rdata"))  ## Plot -- Updated Ash
# load(paste0(datafdr,"dataNHim_withChk_3_sets_PhotoScore0123.rdata"))  ## Individual

load(paste0(datafdr,"dataNHpi_withChk_3Yrs_PhotoScore123_07152021.rdata"))  
load(paste0(datafdr,"dataNHim_withChk_3Yrs_PhotoScore0123_07152021.rdata"))

dataNHpi<-dataNHpi3yrs_C  ##!!!!!

####### These below are the data formatting procedure 
dataNHpi$densityBlades<-ifelse(dataNHpi$densityBlades==0,NA,dataNHpi$densityBlades)  # densityblades as 0 then NA
dataNHpi$popChk <- ifelse(substr(dataNHpi$plotNo, 1, 1) == "Z", substr(dataNHpi$plotNo, 1, 2), "ES")  # Checks VS ES

dataNHpi$withinLoc <- ifelse(as.vector(dataNHpi$femaParLoc) == as.vector(dataNHpi$maleParLoc), 1, 0) # WithinLoc is 1
dataNHpi$development[dataNHpi$development=="#N/A"] <-NA
  str(dataNHpi)

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


### Plot Level----Add blades!!!!
#traits<-c("wetWgtPerM","percDryWgt","dryWgtPerM","densityBlades") #"Ash","AshFDwPM",
# for (col in traits){
#   dataNHpi[,col]<-as.numeric(dataNHpi[,col])
# }

#
  dataNH0<-dataNHpi   ### !!!!! Only the plot traits


#### Individual level to get their experimental factors.
#### The dataNHpi is already filtered for their phenotypic PHOTO SCORE (>2)
  
dataNHim<-dataNHim3yrs_C
dataNHim$line<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "line")
dataNHim$block<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "block")
dataNHim$Year<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "Year")
dataNHim$popChk<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "popChk")
dataNHim$Crosses<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "Crosses")
dataNHim$PhotoScore<-expss::vlookup(dataNHim$plotNo,dict=dataNHpi,lookup_column = "plotNo",result_column = "PhotoScore")
dataNHim<-dataNHim[which(dataNHim$PhotoScore >1),]  # 3969 rows with PhotoScore >1
  str(dataNHim)
  dim(dataNHim)
# 
# traits<-c("bladeLength","bladeMaxWidth","bladeThickness","stipeLength","stipeDiameter")
# for (col in traits){
#   dataNHim[,col]<-as.numeric(dataNHim[,col]) }

#######
####### Averaging the individual measurements to a plot level
library(dplyr)

dataNHim_avg<-matrix(data=0,nrow=372,ncol=1) # 372 is the # of rows in EachTrait
traitInd<-c("bladeLength","bladeMaxWidth","bladeThickness","stipeLength","stipeDiameter")  
   for (i in 1:length(traitInd)){
    dataNHim$Trait<-dataNHim[,traitInd[i]]
    EachTrait<-aggregate(Trait ~ Year+line+block+Crosses, data=dataNHim, FUN=mean, na.action = na.omit)
    colnames(EachTrait)[colnames(EachTrait)=="Trait"]<-trait[i]
    dataNHim_avg<-cbind(dataNHim_avg,EachTrait)
   }

dataNHim_avg<-dataNHim_avg[,c("Year","line","block","Crosses","bladeLength","bladeMaxWidth","bladeThickness","stipeLength","stipeDiameter")]
#dataNHim_avg<-aggregate(cbind(bladeLength,bladeMaxWidth,bladeThickness,stipeLength,stipeDiameter) ~ Year+line+block+Crosses, data=dataNHim, FUN=mean, na.action = na.omit)
### Plot and individual traits
dataNHpi$MergeVar<-paste(dataNHpi$Year,dataNHpi$line,dataNHpi$block,dataNHpi$Crosses,sep="_")
dataNHim_avg$MergeVar<-paste(dataNHim_avg$Year,dataNHim_avg$line,dataNHim_avg$block,dataNHim_avg$Crosses,sep="_")

dataNH1<-merge(dataNHpi,dataNHim_avg,by.x="MergeVar",by.y="MergeVar",all.x=TRUE)
  dataNH1[is.na(as.character(dataNH1$Crosses.y)==as.character(dataNH1$Crosses.y)),]
colnames(dataNH1)[colnames(dataNH1)=="Year.x"]<-"Year"
colnames(dataNH1)[colnames(dataNH1)=="line.x"]<-"line"
colnames(dataNH1)[colnames(dataNH1)=="block.x"]<-"block"
colnames(dataNH1)[colnames(dataNH1)=="Crosses.x"]<-"Crosses"
  colnames(dataNH1)
#######
  

############
#### Between Year Genetic cor for all traits
###asreml.license.activate()
###enter this code CDEA-HECC-CDAH-FIED

  library(asreml)
  
  ### Only the individual traits
  # dataNH0<-dataNHim   ### !!!!!  OR dataNH<-dataNHpi ----> In 5.2 script RUN
  
  dataNH0<-dataNH1  ###!!!!!  Both the plot level and individual level traits are merged now! 
  
  
# ### Yr19, 21
dataNH<-droplevels(dataNH0[dataNH0$Year%in%c(2019,2021),])
scheme<-"Yr19_Yr21"

# dataNH<-droplevels(dataNH0[dataNH0$Year%in%c(2019,2020),])
# scheme<-"Yr19_Yr20"

 # dataNH<-droplevels(dataNH0[dataNH0$Year%in%c(2020,2021),])
 # scheme<-"Yr20_Yr21"


### Checks in the bottom
dataNH_RMchk<-dataNH[!dataNH$crossID=="Check",]
dataNH_chks<-dataNH[dataNH$crossID=="Check",]
dataNH<-rbind(dataNH_RMchk,dataNH_chks)  
  dim(dataNH)
dataNH2<-dataNH

traits<-c("wetWgtPerM","percDryWgt","dryWgtPerM","densityBlades","bladeLength","bladeMaxWidth","bladeThickness","stipeLength","stipeDiameter")
 
modSum<-NULL
cor<-NULL


# Gencor_Yr<-function (data=dataNH,Trait_grm=Trait_grm3){
#   
#   mod4<-asreml(Trait~Year+line+block,
#                random= ~ us(Year):vm(Crosses,Trait_grm3),
#                data = data, maxiter=500, trace=TRUE)
#   #+popChk Yr19_Yr21, specify without popChk in the model
#   return(list(mod=mod4))
# }

Gencor_Yr<-function (data=dataNH,Trait_grm=Trait_grm3){
  
  mod4<-asreml(Trait~Year+line+block+popChk,
               random= ~ us(Year):vm(Crosses,Trait_grm3),
               data = data, maxiter=100, trace=TRUE)
  return(list(mod=mod4))
}



### Yr19_Yr21 is c(1:2,5:length(traits))

for (j in c(1:2,4:length(traits))){
dataNH<-dataNH2
dataNH$Trait<-as.numeric(dataNH[,traits[j]])
dataNH<-droplevels(dataNH[!is.na(dataNH$Trait),])

# if (traits[j]%in%c("dryWgtPerM")){ #,"AshFDwPM" this is for 2019 and 2020 only
#   dataNH$Trait2<-ifelse(dataNH$Year==2019,dataNH$Trait*sqrt(10),dataNH$Trait)  # Yr1 phenotype * sqrt(10)
#   dataNH$Trait<-dataNH$Trait2
# }

Trait_Crosses<-as.character(dataNH$Crosses)
Trait_grm<-All_grm[rownames(All_grm)%in%Trait_Crosses,colnames(All_grm)%in%Trait_Crosses]
    print(dim(dataNH))
    print(sum(is.na(dataNH$Trait)))
    print(dim(Trait_grm))

### Adding the checks into the Trait_grm, all 1s in diagonal,all 0s for others
data<-dataNH    
    #data[data$popChk=="ES",]$Crosses
    droplevels(data[data$popChk=="ES",])$Crosses
    ChkCross<-unique(droplevels(data[!data$popChk=="ES",])$Crosses) # 33 plots of checks, 4 unique ones
    Col0<-matrix(0,nrow=nrow(Trait_grm),ncol=length(ChkCross))
    colnames(Col0)<-ChkCross
    
    Trait_grm2<-cbind(Trait_grm,Col0)
    
    Row0<-matrix(0,nrow=length(ChkCross),ncol=ncol(Trait_grm))
    Chk1<-diag(x=1,nrow=length(ChkCross),ncol=length(ChkCross))
    Row0_Chk1<-cbind(Row0,Chk1)
    rownames(Row0_Chk1)<-ChkCross
    
Trait_grm3<-rbind(Trait_grm2,Row0_Chk1)

## Correct the order Trait_grm3 to match to Crosses
#Trait_grm3<-Trait_grm3[levels(data$Crosses),levels(data$Crosses)]

Trait_grm3<-Trait_grm3[match(levels(data$Crosses),rownames(Trait_grm3)),match(levels(data$Crosses),colnames(Trait_grm3))]
    
  print(nrow(Trait_grm3)==length(unique(droplevels(data$Crosses))))
  identical(rownames(Trait_grm3),levels(data$Crosses))
  
modSum[[j]]<-summary(Gencor_Yr(data=data,Trait_grm=Trait_grm3)$mod)$varcomp
#covarianceA_B/sqrt(varianceA*varianceB)

cor<-c(cor,modSum[[j]][,"component"][2]/sqrt(modSum[[j]][,"component"][1]*modSum[[j]][,"component"][3]))
}

### Q1?????????? Trait_grm3 is 212x212 VS the data/dataNH is 247xntrait
### Add PlotValue 1:nrow(data), Add the NAs to those plotValue
# for (year in c(2020,2021)){
#   data[data$Year==year,]$plotdummy<-1:nrow(data[data$Year==year,])
# }

# mod4<-asreml(Trait~Year+line+block+popChk,
#              random= ~ us(Year):vm(Crosses,Trait_grm3),
#              residual= ~ idh(Year),
#              data = data, maxiter=100, trace=TRUE,na.action=na.method(x="include"))


#:id(units) #:id(plotdummy)
### NOW Assumed error variance to be the same.
### Adding diagnol error variance structure. Within Yr homogeneous/ Between hetero

names(cor)<-traits 
##!! For scheme: Yr19_Yr21, DwpM has non-positive variance str
# names(cor)<-traits[!traits==traits[3]]


# ##!!!!!!!!!!!!!!!!!!!!!!!!
# #### !!! Between Yr, Same Trait Plot, Gen_cor  ### ----> In 5.2 script RUN
# Plotcor<-cor  
# write.csv(cor,paste0("Genetic_Correlations_BetweenYears_PlotLevel_",scheme,"_0712_2021.csv"))
# save(modSum,file=paste0("Genetic_Correlations_BetweenYears_PlotLevel_moddelSummary_",scheme,"_0712_2021.Rdata"))
# 
# ##!!!!!!!!!!!!!!!!!!!!!!!!!
# #### !!! Between Yr, Same Trait Individual, Gen_cor
# Indicor<-cor
# write.csv(cor,paste0("Genetic_Correlations_BetweenYears_IndivLevel_",scheme,"_0712_2021.csv"))
# save(modSum,file=paste0("Genetic_Correlations_BetweenYears_IndivLevel_modelSummary_",scheme,"0712_2021.Rdata"))
# 
# Bothcor<-c(Plotcor,Indicor)
  Bothcor<-cor
write.csv(Bothcor,paste0("Genetic_Correlations_BetweenYears_Plot_Indiv_Level_",scheme,"_0712_2021.csv"))




#############
#############

########### Multi-trait for plot-level +Individual level traits
################################## Adding the checks into Relationship matrix
  ### Transform the traits value by multiply sqrt(10)
dataNH<-dataNHpi
# for (j in 1: length(traits)){
# dataNH$Trait<-dataNH[,traits[j]]
# dataNH<-droplevels(dataNH[!is.na(dataNH$Trait),])
# 
# if (traits[j]%in%c("dryWgtPerM","AshFDwPM")){
#   dataNH$Trait2<-ifelse(dataNH$Year==2019,dataNH$Trait*sqrt(10),dataNH$Trait)  # Yr1 phenotype * sqrt(10)
#   dataNH$Trait<-dataNH$Trait2
#   }
# }

Trait_Crosses<-as.character(dataNH$Crosses)
Trait_grm<-All_grm[rownames(All_grm)%in%Trait_Crosses,colnames(All_grm)%in%Trait_Crosses]
  print(dim(dataNH))
  print(sum(is.na(dataNH$Trait)))
  print(dim(Trait_grm))

### Adding the checks into the Trait_grm, all 1s in diagonal,all 0s for others
data<-dataNH1 ### data has both of the Plot level and Individual level traits

data[data$popChk=="ES",]$Crosses
droplevels(data[data$popChk=="ES",])$Crosses
ChkCross<-unique(droplevels(data[!data$popChk=="ES",])$Crosses) # 33 plots of checks, 4 unique ones
Col0<-matrix(0,nrow=nrow(Trait_grm),ncol=length(ChkCross))
colnames(Col0)<-ChkCross

Trait_grm2<-cbind(Trait_grm,Col0)

Row0<-matrix(0,nrow=length(ChkCross),ncol=ncol(Trait_grm))
Chk1<-diag(x=1,nrow=length(ChkCross),ncol=length(ChkCross))
Row0_Chk1<-cbind(Row0,Chk1)
rownames(Row0_Chk1)<-ChkCross

Trait_grm3<-rbind(Trait_grm2,Row0_Chk1)

print(nrow(Trait_grm3)==length(unique(droplevels(data$Crosses))))

## Store and conver Crosses into numeric form
## Re-name the Crosses numeric! and Each order has to match up with the order in the Gsnp

Trait_grm4<-Trait_grm3 # store out grm3 (grm4 still has actual character as rownames, grm3 will be numeric)
Cross_link<-as.data.frame(rownames(Trait_grm3))
colnames(Cross_link)<-"Cross_chr"

Cross_link$num<-1:nrow(Trait_grm3)
    ChkCross
    
#Cross_link$num<-ifelse(as.character(Cross_link$Cross_chr)%in%ChkCross,999,Cross_link$num0)

rownames(Trait_grm3)<-expss::vlookup(rownames(Trait_grm3),dict=Cross_link,result_column = "num",lookup_column ="Cross_chr" )
colnames(Trait_grm3)<-expss::vlookup(colnames(Trait_grm3),dict=Cross_link,result_column = "num",lookup_column ="Cross_chr" )

data$CrossesNum<-expss::vlookup(data$Crosses,dict=Cross_link,result_column = "num",lookup_column ="Cross_chr" )

data$CrossesNum<-as.factor(data$CrossesNum)
### Order the data$CrossesNum to be matching with the Trait_grm3 
data<-data[order(data$Year,data$CrossesNum),]

### Inverse the Trait_grm3
snpRelMat<-Trait_grm3
Gsnp=solve(snpRelMat+diag(1e-6,length(snpRelMat[,1]),length(snpRelMat[,1])))

# Map elements in the relationship matrix to the phenotypes

nrow(Gsnp)==sum(rownames(Gsnp)==levels(data$CrossesNum))


attr(Gsnp, "INVERSE")=TRUE

asreml.options(extra=30)

#Within plot level traits and Individual traits
modMT1<-asreml(cbind(percDryWgt,dryWgtPerM,AshFDwPM,densityBlades,bladeLength,bladeMaxWidth,bladeThickness,stipeLength,stipeDiameter) ~ trait+trait:line+trait:block+trait:popChk+trait:Year,
             random= ~ us(trait):vm(CrossesNum,Gsnp),
             residual = ~id(units):us(trait),
             data = data, maxiter=100, trace=TRUE)  # Final Used 4 Plot traits +5 Indi traits

# modMT2<-asreml(cbind(dryWgtPerM,densityBlades,AshFDwPM,percDryWgt) ~ trait+trait:line+trait:block+trait:popChk+trait:Year,
#                random= ~ us(trait):vm(CrossesNum,Gsnp),
#                residual = ~id(units):us(trait),
#                data = data, maxiter=100, trace=TRUE)
# 
# modMT3<-asreml(cbind(bladeLength,bladeMaxWidth,bladeThickness,stipeLength,stipeDiameter) ~ trait+trait:line+trait:block+trait:popChk+trait:Year,
#                random= ~ us(trait):vm(CrossesNum,Gsnp),
#                residual = ~id(units):us(trait),
#                data = data, maxiter=100, trace=TRUE)
# 
# modMT4<-asreml(cbind(percDryWgt,dryWgtPerM,Ash,AshFDwPM,densityBlades,bladeLength,bladeMaxWidth,bladeThickness) ~ trait+trait:line+trait:block+trait:popChk+trait:Year,
#                random= ~ us(trait):vm(CrossesNum,Gsnp),
#                residual = ~id(units):us(trait),
#                data = data, maxiter=100, trace=TRUE)
### Use id(CrossesNum) for 4 plot traits runs
### Use vm(CrossesNum,Gsnp) for 2 plot traits dryWgtPerM and densityBlades runs fine
### Use vm(CrossesNum,Gsnp) for 3 plot traits (+AshFDWpM), warnings: convergence "modify to make positive)
### Use vm(CrossesNum,Gsnp)  for 4 plot traits (+AshFDWpM,wetWgtPerM) gives "Log-likelihood not converged"
### Use vam(CrossesNum,Gsnp) for 4 plot traits (+AshFDWpM,percDryWgt) 
### Use vm(CrossesNum,Gsnp) for 2 plot traits+ Individual traits, runs

# 
# mod3<-asreml(wetWgtPerM~Year+line+block+popChk,
#              random= ~ us(Year):vm(CrossesNum,Trait_grm3),
#              data = data, maxiter=100, trace=TRUE)


#### Getting the genetic cor
GenVar1<-summary(modMT1)$varcomp
write.csv(GenVar1,"Plot4traits_Indivial_MultiTrait_Genetic_var_05272021.csv") #and_

  colnames(GenVar1)
## Extract rows with extra parts

Gen_varcov<-GenVar1[grepl("trait:vm",rownames(GenVar1)),]
## RM names with extra parts
rownames(Gen_varcov)<-stringr::str_split_fixed(rownames(Gen_varcov),"_",2)[,2]

#Mtraits<-c("wetWgtPlot","dryWgtPerM","AshFDwPM","densityBlades")

### The varcomp is in long format, convert to wide format
Gen_varcov$var1<-stringr::str_split_fixed(rownames(Gen_varcov),":",2)[,1]
Gen_varcov$var2<-stringr::str_split_fixed(rownames(Gen_varcov),":",2)[,2]
Gen_varcov$cor<-Gen_varcov$component
Gen_varcov<-Gen_varcov[order(rownames(Gen_varcov)),]
Gen_varcov_wide<-tidyr::spread(Gen_varcov[,c("var1","var2","component")],var2,component)

write.csv(Gen_varcov_wide,"Plot4traits_Individual_Gen_varcov_wide_05272021.csv")

### Reorganize: Fill in NAs
rownames(Gen_varcov_wide)<-Gen_varcov_wide$var1
Gen_varcov_wide<-Gen_varcov_wide[,-1]
Gen_varcov_imp<-matrix(nrow=nrow(Gen_varcov_wide),ncol=ncol(Gen_varcov_wide))
rownames(Gen_varcov_imp)<-rownames(Gen_varcov_wide)
colnames(Gen_varcov_imp)<-colnames(Gen_varcov_wide)

for (i in 1:nrow(Gen_varcov_wide)){
  for (j in 1:ncol(Gen_varcov_wide)){
    Gen_varcov_imp[i,j]<-ifelse(is.na(Gen_varcov_wide[i,j]),Gen_varcov_wide[j,i],Gen_varcov_wide[i,j])
  }
}

isSymmetric(Gen_varcov_imp)
Gen.cor<-cov2cor(Gen_varcov_imp)
orders<-c("percDryWgt","dryWgtPerM","AshFDwPM","densityBlades","bladeLength","bladeMaxWidth","bladeThickness","stipeLength","stipeDiameter")

Gen.cor<-Gen.cor[match(orders,rownames(Gen.cor)),match(orders,colnames(Gen.cor))]

write.csv(Gen.cor,"Plot4traits_Individual_5traits_Multi_trait_Genetic_cor_05272021.csv") #

### If need to Make varcov symmetric
# library(gdata)
# lowerTriangle(Gen_varcov_wide)<-upperTriangle(Gen_varcov_wide,byrow=TRUE)
# Gen_varcov_wide<-as.matrix(Gen_varcov_wide)  







#######COMPARE?? Multi-Trait using the dBLUPs!!
source("/local/workdir/mh865/GCA_SCA/OneTime1920/code/is.positive.definite.R")
source("/local/workdir/mh865/GCA_SCA/OneTime1920/code/is.square.matrix.R")
source("/local/workdir/mh865/GCA_SCA/OneTime1920/code/is.symmetric.matrix.R")
# 

# Y<-dataNHpi

#### Use the de-regressed BLUPs for multi-trait analysis?
datafdr<-paste0(WD,"OneTime1920/data/")
load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_Individuals_level_WithinYear.Rdata")) ##!!!

Y<-WithinYr_Both_dBLUPs[,!colnames(WithinYr_Both_dBLUPs)%in%c("Row.names","Year.y","plotNo.y","Crosses.y")]
rownames(Y)<-Y$plotNo.x
Y$Crosses<-as.factor(Y$Crosses.x)
Y$Year<-Y$Year.x

is.positive.definite(outCovComb4_dipOrder)
outCovComb4_dipOrder[1:4,1:5]

Amat<-outCovComb4_dipOrder[rownames(outCovComb4_dipOrder)%in%as.character(Y$Crosses),colnames(outCovComb4_dipOrder)%in%as.character(Y$Crosses)]

snpRelMat<-Amat
Gsnp=solve(snpRelMat+diag(1e-6,length(snpRelMat[,1]),length(snpRelMat[,1])))

# Map elements in the relationship matrix to the phenotypes
  nrow(Gsnp)==length(unique(Y$Crosses))
  sum(rownames(Gsnp)==levels(Y$Crosses))
  sum(colnames(Gsnp)==levels(Y$Crosses))

Gsnp2<-Gsnp[match(levels(Y$Crosses),rownames(Gsnp)),match(levels(Y$Crosses),colnames(Gsnp))]
Gsnp<-Gsnp2



attr(Gsnp, "INVERSE")=TRUE

for (t in c("dryWgtPerM","AshFDwPM")){
  Y$Trait<-ifelse(Y$Year==2019,Y[,t]*sqrt(10),Y[,t])  # Yr1 phenotype * sqrt(10)
  colnames(Y)[colnames(Y)=="Trait"]<-paste0(t,"_sqrt10")
}


#### Multi-Trait that worked, for morphology traits dBLUPs
modMTM <- asreml(cbind(bladeLength,bladeMaxWidth,bladeThickness,stipeLength,stipeDiameter) ~ trait, 
                 random= ~ us(trait):vm(Crosses,Gsnp),
                 residual = ~id(units):us(trait),
                 data = Y, maxiter=50, trace=TRUE,extra=10) #

# Did not work for dBLUPs of plot level traits
# modMTM2 <- asreml(cbind(dryWgtPerM_sqrt10,wetWgtPerM) ~ trait, 
#                  random= ~ us(trait):vm(Crosses,Gsnp),
#                  residual = ~id(units):us(trait),
#                  data = Y, maxiter=50, trace=TRUE)

#### May not enough E variance from dBLUPs
#### Trouble shoot: Run uni-variate on dBLUPs (if it converges)


#generating estimates of Crosses performance in each trait
predMTM = predict(modMTM, classify = "trait:Crosses", trace=F)

# pulling out the predictions and comparing to the true simulated values
pMTM=predMTM$pvals$predicted.value


#### Getting the Genetic cor values

GenVar<-summary(modMTM)$varcomp
write.csv(GenVar,"Individual_MultiTrait_Genetic_var.csv")

getwd()
##Manually replaced extra char strings !!!
Indi_Genetic_Cor<-read.csv(paste0(datafdr,"Individual_MultiTrait_Genetic_var.csv"),row.names=1)
head(Indi_Genetic_Cor)

grepl(":",Indi_Genetic_Cor$X.1)
Indi_traits<-c("bladeLength","bladeMaxWidth","bladeThickness","stipeLength","stipeDiameter")
rownames(Indi_Genetic_Cor)[1:(3*length(Indi_traits))]
  #coVarA_B/sqrt(VarA*VarB)

for (i in 1:length(Indi_traits)){
  trait<-Indi_traits[i]
  traitname<-Indi_Genetic_Cor$X.1
  Itself[i]<-Indi_Genetic_Cor[,"component"][grepl(paste0(trait,":",trait),traitname)]
  WithOther[i]<-Indi_Genetic_Cor[,"component"][grepl(paste0(trait,":",))]
}