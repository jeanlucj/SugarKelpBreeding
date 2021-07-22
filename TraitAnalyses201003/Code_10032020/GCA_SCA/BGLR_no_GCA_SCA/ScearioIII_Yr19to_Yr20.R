#ScenarioIII_3.d_Yr19To20.R

#ScenarioIII_3.c_Yr19To21.R

### Yr2019 predict Yr2020

rm(list=ls())
#
WD<-"/Users/maohuang/Desktop/Kelp/Simulation_Study/SugarKelpBreeding/TraitAnalyses201003/"
#WD<-"/local/workdir/mh865/GCA_SCA/OneTime192021/"  # run in terminal

datafdr<-paste0(WD,"data/")

####1. To get the experimental factors column Info
load(paste0(datafdr,"dataNHpi_withChk_3Yrs_PhotoScore123_07152021.rdata"))   ## Plot
dataNHpi<-dataNHpi3yrs_C

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


Y1<-dataNHpi

colKeep<-c("crossID","Crosses","femaPar","femaParLoc","malePar","maleParLoc","plotNo","Region","popChk","line","block","Year","PhotoScore","dryWgtPerM") #,"AshFreedryWgtPerM"
Y2<-droplevels(Y1[,colKeep])
  head(Y2)
Y<-Y2 # Both Years
Y<-droplevels(Y[Y$popChk=="ES",])
ExpCol<-Y[,colnames(Y)%in%c("crossID","Crosses","femalPar","femaParLoc","malePar","maleParLoc","plotNo","Region","popChk","line","block","Year","PhotoScore")]

##### !!!
pop<-"Yr19to20"  ##!!!
TPYr<-2019      # TPYr   
NAYr<-2020 
#####!!!

#### 2. Upload the dBLUPs (estimated WithinYr 2019 and 2020, respectively)
load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_Individuals_level_WithinYear_AddBD.Rdata"))  #2019 and 2020 

load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_n_Indiv_WithinYr2021_0715_2021.Rdata")) # 2021

### Best to merge the Yr2021 data to the other two years data
### This is just to make the Yr2021 colnames the same as the other two years

WithinYr_Both_dBLUPs<-WithinYr_Both_dBLUPs[,!colnames(WithinYr_Both_dBLUPs)%in%c("Ash","AshFDwPM")]
colnames(WithinYr21_dBLUPs)[!colnames(WithinYr21_dBLUPs)%in%colnames(WithinYr_Both_dBLUPs)]
colnames(WithinYr_Both_dBLUPs)[!colnames(WithinYr_Both_dBLUPs)%in%colnames(WithinYr21_dBLUPs)]

WithinYr_Both_dBLUPs<-WithinYr_Both_dBLUPs[,colnames(WithinYr21_dBLUPs)]
identical(colnames(WithinYr_Both_dBLUPs),colnames(WithinYr21_dBLUPs))

WithinYr_Both_dBLUPs<-rbind(WithinYr_Both_dBLUPs,WithinYr21_dBLUPs)

rownames(WithinYr_Both_dBLUPs)<-WithinYr_Both_dBLUPs$Row.names
WithinYr_Both_dBLUPs<-WithinYr_Both_dBLUPs[,!colnames(WithinYr_Both_dBLUPs)=="Row.names"]

WithinYr_Both_dBLUPs<-droplevels(WithinYr_Both_dBLUPs)
  dim(WithinYr_Both_dBLUPs) 
  head(WithinYr_Both_dBLUPs)
identical(as.character(rownames(WithinYr_Both_dBLUPs)),as.character(WithinYr_Both_dBLUPs$plotNo.x))

traits<-colnames(WithinYr_Both_dBLUPs)[!colnames(WithinYr_Both_dBLUPs)%in%c("Row.names","Crosses.x","Crosses.y","plotNo.x","plotNo.y","Year.x","Year.y")]
#2021 does not have Ash; Its percDryWgt is still 0s

traits<-traits[!traits%in%c("Ash","AshFDwPM")]
  print(traits)
save(WithinYr_Both_dBLUPs,file=paste0(datafdr,"dBLUPs_3_yrs_cal_WithinYr_",pop,".Rdata"))

##### 

### 3. Merge the experimental factors col to drBLUP that is used for Training Y1.0
Y1.0<-droplevels(ExpCol[ExpCol$Year%in%TPYr,])  # TPYr 2019
CrossBLUEYr<-WithinYr_Both_dBLUPs[WithinYr_Both_dBLUPs$Year.x%in%TPYr,]
CrossBLUEYr<-droplevels(CrossBLUEYr[,!colnames(CrossBLUEYr)%in%c("Year.x","Year.y","Crosses.y","plotNo.x","plotNo.y")])

Y1<-merge(CrossBLUEYr,Y1.0,by.x="row.names",by.y="plotNo",all.x=TRUE)  ### If by "Crosses', there will be dup in 2019,2020
  dim(Y1)
  head(Y1)
Y1$plotNo<-Y1$Row.names
Y1<-Y1[,!colnames(Y1)=="Row.names"]

#### This is the prediction PP
Y2.0<-ExpCol[ExpCol$Year==NAYr,]  # NAYr 2021
CrossBLUEYr<-WithinYr_Both_dBLUPs[WithinYr_Both_dBLUPs$Year.x==NAYr,]
CrossBLUEYr<-CrossBLUEYr[,!colnames(CrossBLUEYr)%in%c("Year.x","Year.y","Crosses.y","plotNo.x","plotNo.y")]
Y2<-merge(Y2.0,CrossBLUEYr,all.x=TRUE,by.x="Crosses",by.y="Crosses.x")
  dim(Y2)
### Y1 has one extra col
Y1<-Y1[,colnames(Y1)%in%colnames(Y2)]
colnames(Y2)%in%colnames(Y1)
colnames(Y1)%in%colnames(Y2)
Y1<-Y1[,colnames(Y2)]  
  identical(colnames(Y1),colnames(Y2))
Yrbind<-rbind(Y1,Y2)
Ydata<-Y   # Save out the original Y, for in case use
Y<-Yrbind  # Now Y is updated with its BLUEs from within 2019, within 2020  and also WithinYr21
  dim(Y)
  str(Y)

# 
# # # ##########DO it ONCE !!!!
# # # ### Load the A matrix, outCovComb
load(paste0(datafdr,"outCovComb4_Mix_Conden_0712_2021.Rdata")) ##!!!

A_950<-outCovComb4_MixOrder
IDs<-as.character(unique(Y$Crosses))  ### RM checks??
A<-as.matrix(A_950[IDs,IDs])

####### Some plots had data, but are the same duplicated cross between years.
####### Need to expand the A matrix to make it fit all those 250 crosses (Y file (6 plots contained duplicated Crosses))
####### This is just to expand the A and match its order to Y!!!
Y$plotNo<-as.factor(Y$plotNo)
Ydup<-as.character(droplevels(Y$Crosses[duplicated(Y$Crosses)])) # Crosses in Y duplicated
Ydup
Y[Y$Crosses%in%Ydup,1:9]

YdupPlotNo<-as.character(droplevels(Y[which(duplicated(Y$Crosses)),]$plotNo)) ### plotNo for these dups
  str(YdupPlotNo)
YdupCrosses<-as.character(droplevels(Y[which(Y$plotNo%in%YdupPlotNo),]$Crosses))
  str(YdupCrosses)

YdupPlotCross<-cbind(YdupPlotNo,YdupCrosses)  

### The Dup within YdupPlotCross is due to dup within a year 2021
YRMdupPlot<-YdupPlotCross[duplicated(YdupPlotCross[,2]),1] #first Col is the YdupPlotNo

## WithinYr DupCrosses will be the same dBLUPs value, so remove one set
Y<-Y.y<-droplevels(Y[!Y$plotNo%in%YRMdupPlot,])   # RM the dup plot within a year (2021), because their dBLUP is the same
## Update the list of dupplot and dupCross

YdupCrossesUpdate<-YdupPlotCross[!YdupPlotCross[,1]%in%YRMdupPlot,2] # 2nd Col is the YdupCross
# YdupPlotNo
YdupPlotNoUpdate<-YdupPlotCross[!YdupPlotCross[,1]%in%YRMdupPlot,1]
  str(YdupCrosses)
  str(YdupCrossesUpdate)
rownames(A)[!rownames(A)%in%as.character(Y.y$Crosses)]  #2020 has 2 dup crosses

A.DupxDup<-A[rownames(A)%in%YdupCrossesUpdate,colnames(A)%in%YdupCrossesUpdate]
A.DupxAll<-A[rownames(A)%in%YdupCrossesUpdate,]
A.AllxDup<-A[,colnames(A)%in%YdupCrossesUpdate]
  identical(rownames(A.DupxDup),rownames(A.DupxAll))
  identical(rownames(A.DupxAll),colnames(A.AllxDup))

Anew.part1<-cbind(A,A.AllxDup)
Anew.part2<-cbind(A.DupxAll,A.DupxDup)
  identical(rownames(A),rownames(A.AllxDup))
  identical(colnames(A),colnames(A.DupxAll))
  identical(colnames(Anew.part1),colnames(Anew.part2))

Anew<-rbind(Anew.part1,Anew.part2)
  dim(Anew)

# The old A part, look up plotNo in the Y data file that has removed the YRMdupPlotUpdate
Y.y<-Y[!Y$plotNo%in%YdupPlotNoUpdate,]
Y.y$plotNo<-as.character(Y.y$plotNo) #Note cannot use Y.y directly in model as this becomes character
  str(Y.y)
ToPlotNo1<-expss::vlookup(lookup_value=rownames(Anew.part1),dict=droplevels(Y.y),result_column="plotNo",lookup_column="Crosses")
  str(ToPlotNo2)
ToPlotNo2<-YdupPlotNoUpdate  # PlotNo for these duplicates
# ToPlotNo2<-expss::vlookup(lookup_value=rownames(Anew.part2),dict=droplevels(Y2),result_column="plotNo",lookup_column="Crosses")
#ToPlotNo<-expss::vlookup(lookup_value=rownames(Anew),dict=Y,result_column="plotNo",lookup_column="Crosses")
ToPlotNo<-c(as.character(ToPlotNo1),ToPlotNo2)
  str(ToPlotNo)
rownames(Anew)<-ToPlotNo
CrossesInAnew<-colnames(Anew)
colnames(Anew)<-rownames(Anew)

Anew<-Anew[match(as.character(Y$plotNo),rownames(Anew)),match(as.character(Y$plotNo),colnames(Anew))]

#Y$plotNo<-as.factor(Y$plotNo)
identical(rownames(Anew),as.character(Y$plotNo))
identical(colnames(Anew),as.character(Y$plotNo))
dim(Anew)
str(Y$plotNo)
length(CrossesInAnew)
pop
save(A,CrossesInAnew,Anew,Y,file=paste0(datafdr,"A_n_Y_",pop,"_allES_matchedtoY.Rdata"))

# # ########## DO it ONCE



#### Load Anew that's already matched order to Y(250 ES indi)
#load("/local/workdir/mh865/GCA_SCA/GBLUPOnly_NoGCA_SCA/A_250Indi_allES_matchedtoY.Rdata")
WD<-"/local/workdir/mh865/GCA_SCA/OneTime192021/"  # run in terminal
datafdr<-paste0(WD,"data/")

### !!!
pop<-"Yr19to20"
TPYr<-2019      # TPYr   
NAYr<-2020 

load(paste0(datafdr,"dBLUPs_3_yrs_cal_WithinYr_",pop,".Rdata"))
traits<-colnames(WithinYr_Both_dBLUPs)[!colnames(WithinYr_Both_dBLUPs)%in%c("Year.x","plotNo.x","Crosses.x","Year.y","plotNo.y","Crosses.y")]
print(traits)

load(paste0(datafdr,"A_n_Y_",pop,"_allES_matchedtoY.Rdata"))

ETA<-list(list(K=Anew,model="RKHS"))
reps<-20 # !!!
ntraits<-length(traits)  # !!!
head(Y)  
str(Y)
identical(rownames(Anew),as.character(Y$plotNo)) 

cor<-matrix(nrow=reps,ncol=length(traits))
colnames(cor)<-traits

for (j in 1:length(traits)){
  Coltrait<-traits[j]
  
  for (i in 1:reps){
    setwd(paste0(WD,pop,"_output/"))
    
    dir.create(paste0(Coltrait,"_ydrBLUP_noLocRep",i))
    savepath<-paste0(pop,"_output/",Coltrait,"_ydrBLUP_noLocRep",i,"/")
    
    
    Y$BLUE_Trait<-Y[,colnames(Y)==Coltrait]
    
    y<-yBLUE<-Y[,"BLUE_Trait"]   # phenotypes column  !!!
    
    yNA<-y
    
    testing<-which(Y$Year==NAYr)  ## PP year
    yNA[testing]<-NA
    
    fm<-BGLR::BGLR(y=yNA,
                   ETA=ETA,
                   nIter=80000,
                   burnIn=60000,
                   saveAt=paste0(WD,savepath,pop,"_rep",i),
                   verbose=TRUE)
    save(fm,file=paste0(WD,savepath,"fm_",pop,"_rep",i,".rda"))
    
    yPred<-fm$ETA[[1]]$u  # equals to fm$yHat
    predict<-data.frame(testing,
                        Crosses=Y$Crosses[testing],
                        yBLUE=yBLUE[testing],
                        yPred=yPred[testing],
                        yHat=fm$yHat[testing])
    predict<-droplevels(predict)
    cor[i,j]<-cor(predict$yBLUE,predict$yPred,use="complete")
    write.table(predict,file=paste0(WD,savepath,Coltrait,"_predictions_rep",i,".csv")) 
  }
  rm(fm) 
  unlink("*.dat")
  print(colMeans(cor))  
}

standard_err<-function(x){sd(x,na.rm=TRUE)/sqrt(length(na.omit(x)))}
stderr<-apply(cor,2,standard_err)
cormean<-colMeans(cor)
cor_std<-rbind(cormean,stderr)
rownames(cor_std)<-c("corMean","StdErr")
colnames(cor_std)<-traits

write.csv(cor,paste0(WD,pop,"_output/","cor_",pop,"_",length(traits),"_traits_","_using_ES_dBLUP_0715_72021.csv"))
write.csv(cor_std,paste0(WD,pop,"_output/","cor_",pop,"_",length(traits),"_traits_","_using_ES_dBLUP_Mean_0715_2021.csv"))


