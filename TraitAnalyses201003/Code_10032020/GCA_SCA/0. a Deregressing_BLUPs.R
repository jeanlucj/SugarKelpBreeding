###############################################  
#### Making dBLUPs within Year and between Years
rm(list=ls())
WD<-"/Users/maohuang/Desktop/Kelp/GCA_SCA/" # local
datafdr<-"/Users/maohuang/Desktop/Kelp/Simulation_Study/SugarKelpBreeding/TraitAnalyses201003/data/"

load(paste0(datafdr,"dataNHpi_withChk_3Yrs_PhotoScore123_07152021.rdata"))  
load(paste0(datafdr,"dataNHim_withChk_3Yrs_PhotoScore0123_07152021.rdata"))

dataNHpi<-dataNHpi3yrs_C  ##!!!!!
#dataNHpi<-dataNHpi3yrs_C[!dataNHpi3yrs_C$PhotoScore<2,]  # 10 plots less for 2021

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

###!!!! Plot level
dataNHpi<-droplevels(dataNHpi)
traits<-c("wetWgtPerM", "percDryWgt","dryWgtPerM","densityBlades")
Plotdata<-TRUE   ### !!!!!
dataNH<-dataNHpi   ### !!!!!

####!!!! Individual level to get their experimental factors.
#### The dataNHpi is already filtered for their phenotypic PHOTO SCORE (>2)
dataNHim3yrs_C$Year<-as.factor(dataNHim3yrs_C$Year) ### !!!!

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

traits<-c("bladeLength","bladeMaxWidth","bladeThickness","stipeLength","stipeDiameter")
dataNH<-dataNHim   ### !!!!!
Plotdata<-FALSE  # means Individual data
#######!!!!
#######


#######
library(sommer)
  dim(dataNH)
  
  
###################
### Both Years !!!!
WithinYear<-FALSE
Yr<-c(2019,2020)
dataNH<-droplevels(dataNH[dataNH$Year%in%Yr,])
dataNH<-dataNH[order(dataNH$plotNo),]
scheme<-"Both"


### Both Years !!!!
WithinYear<-FALSE
Yr<-c(2020,2021)
dataNH<-droplevels(dataNH[dataNH$Year%in%Yr,])
dataNH<-dataNH[order(dataNH$plotNo),]
scheme<-"Yr20_21"


###################
### Within Year
  # 2020_273 is the same as 2020_219
 # comX<-c("SL18-OI-15-FG1xSA18-CB-2-MG3","SL18-LD-13-FG2xSL18-LD-13-MG2","SL18-LD-2-FG2xSL18-SF-19-MG1","SL18-NL-2-FG3xSA18-CB-4-MG1","SL18-NL-3-FG1xSA18-CB-7-MG2")
#
Yr<-2021  ###!!!!   ### THIS INCLUDED PhotoScore 1 plots
scheme<-"2021"

#
Yr<-2020  ###!!!!
scheme<-"2020"
#
Yr<-2019  ###!!!!
scheme<-"2019"

WithinYear=TRUE
dataNH<-dataNH[dataNH$Year==Yr,]
dataNH<-dataNH[order(dataNH$plotNo),]
##############


  dim(dataNH)
  head(dataNH)
  
  if(scheme=="Yr20_21" & Plotdata==TRUE){
    traits<-traits[c(1,3:length(traits))]
  }else{
    traits<-traits
  }  
  
if (WithinYear==FALSE){
  ### Both Years
  NumCrosses<-length(unique(dataNH$Crosses))
  dBLUPs<-BLUPs<-matrix(nrow=NumCrosses,ncol=length(traits))
  H2<-NULL
  convInfo<-vector()
 
  # if(scheme=="Yr20_21" & Plotdata==TRUE){
  #   options<-c(1,3:length(traits))
  # }else{
  #   options<-1:length(traits)
  # }
   
for (j in 1:length(traits)){
  Coltrait<-traits[j]
  dataNH$Trait<-dataNH[,Coltrait]  ### !!!!! TraitBLUE
  data<-dataNH
  data<-droplevels(data)
  
  tmp.mod <- mmer(fixed = Trait ~ popChk+Year, 
                  random = ~block + line + block:line+Crosses, 
                  data = data)
  
  PEV <- diag(tmp.mod$PevU$`Crosses`$Trait)
  varG <- as.numeric(tmp.mod$sigma$`Crosses`)
  tmp.BLUPs <- tmp.mod$U$`Crosses`$Trait  ##blups
  names(tmp.BLUPs)<-stringr::str_replace(names(tmp.BLUPs),"Crosses","")
  tmp.dBLUPs<-tmp.BLUPs / (1 - PEV/varG)  ##deregressed blups checks will have "inf" values
  print(identical(names(tmp.BLUPs),names(tmp.dBLUPs)))
  CheckCrosses<-unique(droplevels(data[!data$popChk=="ES",]$Crosses))  ## 4 Crosses
  CheckCrosses
  tmp.BLUPs[names(tmp.BLUPs)%in%CheckCrosses]<-NA
  tmp.dBLUPs[names(tmp.dBLUPs)%in%CheckCrosses]<-NA
  
  BLUPs[,j] <- tmp.BLUPs  
  dBLUPs[,j] <- tmp.dBLUPs  
  
  #H2 <- rbind(H2, pin(tmp.mod, h2 ~ V2 / ( V2 + V4)))
  H2<-rbind(H2,h2.fun(tmp.mod,data=data,gTerm="Crosses"))
  convInfo[j] <- tmp.mod$convergence
  }
  
} else if(WithinYear==TRUE){
  ### Within Year
  NumCrosses<-length(unique(dataNH$Crosses))
  dBLUPs<-BLUPs<-matrix(nrow=NumCrosses,ncol=length(traits))
  H2<-NULL
  convInfo<-vector()

  for (j in 1:length(traits)){
    Coltrait<-traits[j]
    dataNH$Trait<-dataNH[,Coltrait]  ### !!!!! TraitBLUE
    data<-droplevels(dataNH)

   
    ###2019, 2020 model
    if(Yr%in%c(2019,2020)){
      tmp.mod <- mmer(fixed=Trait~ popChk, 
                      random = ~block + line + block:line+Crosses, 
                      data = data)
    }else if (Yr==2021){
      tmp.mod <- mmer(Trait~1,
                      random = ~ line+block+line:block+Crosses,
                      data = data)
    }    ###2021 model

    PEV <- diag(tmp.mod$PevU$`Crosses`$Trait)
    varG <- as.numeric(tmp.mod$sigma$`Crosses`)
    tmp.BLUPs <- tmp.mod$U$`Crosses`$Trait  ##blups
    names(tmp.BLUPs)<-stringr::str_replace(names(tmp.BLUPs),"Crosses","")
    tmp.dBLUPs<-tmp.BLUPs / (1 - PEV/varG)  ##deregressed blups checks will have "inf" values
    print(identical(names(tmp.BLUPs),names(tmp.dBLUPs)))
    
    ### 2021 No chks
    if(Yr%in%c(2019,2020)){
    CheckCrosses<-unique(droplevels(data[!data$popChk=="ES",]$Crosses))  ## 4 Crosses
    CheckCrosses
    tmp.BLUPs[names(tmp.BLUPs)%in%CheckCrosses]<-NA
    tmp.dBLUPs[names(tmp.dBLUPs)%in%CheckCrosses]<-NA
    } else if (Yr==2021){
      tmp.BLUPs<-tmp.BLUPs
    }
    
    BLUPs[,j] <- tmp.BLUPs  
    dBLUPs[,j] <- tmp.dBLUPs  
    
      ##OR H2 <- rbind(H2, pin(tmp.mod, h2 ~ V2 / ( V2 + V4)))
    H2<-rbind(H2,h2.fun(tmp.mod,data=data,gTerm="Crosses"))
    convInfo[j] <- tmp.mod$convergence
  }
  
}

# ########## Trouble shooting??  
# ### 2021 data for percDryWgt has no line effects. BLUPs are 0??????????
#   Coltrait<-traits[j]
#   dataNH$Trait<-dataNH[,Coltrait]  ### !!!!! TraitBLUE
#   data<-dataNH
#   data<-droplevels(data)
#   fit2021<- lme4::lmer(Trait ~ line*block + (1|Crosses), data=data)
# ################################

rownames(H2)<-traits 
  H2
  convInfo

  if(Plotdata==TRUE){
    save(H2,convInfo,file=paste0(datafdr,scheme,"_PlotdeRegress_sommer_h2.fuc_convergence.Rdata")) #/Plots_
  }else if(Plotdata==FALSE){
    save(H2,convInfo,file=paste0(datafdr,scheme,"_IndideRegress_sommer_h2.fuc_convergence.Rdata")) #/Indi_
  }
 
colnames(BLUPs)<-colnames(dBLUPs)<-traits
rownames(BLUPs)<-rownames(dBLUPs)<-names(tmp.dBLUPs)
  head(BLUPs)
  head(dBLUPs)
  dim(BLUPs)
  dim(dBLUPs)

if(WithinYear==TRUE){
  if (Yr %in% c(2019,2020)){
    BLUPs<-BLUPs[!rownames(BLUPs)%in%CheckCrosses,]     #RM the checks crosses
    dBLUPs<-dBLUPs[!rownames(dBLUPs)%in%CheckCrosses,]
  }else if (Yr==2021){
    BLUPs<-BLUPs
    dBLUPs<-dBLUPs
  }
}else if (WithinYear==FALSE){
  BLUPs<-BLUPs[!rownames(BLUPs)%in%CheckCrosses,]     #RM the checks crosses
  dBLUPs<-dBLUPs[!rownames(dBLUPs)%in%CheckCrosses,]
}  
  
####Yr2021/Yr2019/YrBoth are all named "Yr21_dBLUPs"
Yr21_dBLUPs<-as.data.frame(dBLUPs)  ##!!!! rownames are Crosses

# dataNHpi$Crosses<-as.character(dataNHpi$Crosses) 
# dataNHpi$Year<-as.numeric(dataNHpi$Year)
#dataNH is already the subset for that WithinYear scheme Yr
if (WithinYear==TRUE){
  Yr21_dBLUPs$Year<-Yr
  Yr21_dBLUPs$plotNo<-expss::vlookup(rownames(Yr21_dBLUPs),dict=dataNH,lookup_column = "Crosses",result_column = "plotNo")
  Yr21_dBLUPs$Crosses<-rownames(Yr21_dBLUPs)
  rownames(Yr21_dBLUPs)<-Yr21_dBLUPs$plotNo
  
}else if (WithinYear==FALSE){
  Yr21_dBLUPs<-Yr21_dBLUPs
}

  dim(Yr21_dBLUPs)
  head(Yr21_dBLUPs)

 head(Yr21_dBLUPs)

    if (Plotdata==TRUE){
      save(Yr21_dBLUPs,file=paste0(datafdr,"Deregressed_BLUPs_ESplots_plotlevel_WithinYr",scheme,"_0715_2021.Rdata"))
   }else if(Plotdata==FALSE){
      save(Yr21_dBLUPs,file=paste0(datafdr,"Deregressed_BLUPs_ESplots_Indivlevel_WithinYr",scheme,"_07152021.Rdata"))
    }

 
##!!! Combine the plot level and individual level together 
  load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plotlevel_WithinYr",scheme,"_0715_2021.Rdata"))
WithinYr_Plot_dBLUPs<-Yr21_dBLUPs

  load(paste0(datafdr,"Deregressed_BLUPs_ESplots_Indivlevel_WithinYr",scheme,"_07152021.Rdata"))
  #rownames(Yr21_dBLUPs)<-Yr21_dBLUPs$Crosses
  #Yr21_dBLUPs<-Yr21_dBLUPs[,!colnames(Yr21_dBLUPs)%in%c("Year","plotNo","Crosses")]
WithinYr_Indi_dBLUPs<-Yr21_dBLUPs
  

if (WithinYear==TRUE){
WithinYr21_dBLUPs<-merge(WithinYr_Plot_dBLUPs,WithinYr_Indi_dBLUPs,by.x="row.names",by.y="row.names",all.x=TRUE)
}else if (WithinYear==FALSE){
WithinYr21_dBLUPs<-merge(WithinYr_Plot_dBLUPs,WithinYr_Indi_dBLUPs,by.x="row.names",by.y="row.names",all.x=TRUE) 
}
  head(WithinYr21_dBLUPs)
  dim(WithinYr21_dBLUPs)
save(WithinYr21_dBLUPs,file=paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_n_Indiv_WithinYr",scheme,"_0715_2021.Rdata"))

##### This is to combine Yr19 and Yr20, and save things out as their original names
scheme<-"2019"
load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_n_Indiv_WithinYr",scheme,"_0715_2021.Rdata"))
WithinYr_19_dBLUPs<-WithinYr21_dBLUPs

Yr<-"2020"
load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_n_Indiv_WithinYr",scheme,"_0715_2021.Rdata"))
WithinYr_20_dBLUPs<-WithinYr21_dBLUPs
  identical(colnames(WithinYr_20_dBLUPs),colnames(WithinYr_19_dBLUPs))

WithinYr_Both_dBLUPs<-rbind(WithinYr_19_dBLUPs,WithinYr_20_dBLUPs)
save(WithinYr_Both_dBLUPs,file=paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_Individuals_level_WithinYear_AddBD.Rdata"))
  rm(list="WithinYr21_dBLUPs")


# This is the dBLUPs that's done using both yr2019 and yr2020 data, just to resave to match old names  
scheme<-"Both"
load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_n_Indiv_WithinYr",scheme,"_0715_2021.Rdata"))
Both_dBLUPs<-WithinYr21_dBLUPs
  rm(list="WithinYr21_dBLUPs")
save(Both_dBLUPs,file=paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_Individuals_level_overTwoYears_AddBD.Rdata"))





# ####
# Yr19_dBLUPs<-as.data.frame(dBLUPs)  ##!!!!
# Yr19_dBLUPs$Year<-Yr
#   dim(Yr19_dBLUPs)
# Yr19_dBLUPs$plotNo<-expss::vlookup(rownames(Yr19_dBLUPs),dict=dataNHpi[dataNHpi$Year==Yr,],lookup_column = "Crosses",result_column = "plotNo")
# Yr19_dBLUPs$Crosses<-rownames(Yr19_dBLUPs)
# rownames(Yr19_dBLUPs)<-Yr19_dBLUPs$plotNo
#   Yr19_dBLUPs[Yr19_dBLUPs$Crosses=="SL18-LD-13-FG2xSL18-LD-13-MG2",]
# ###
# Yr20_dBLUPs<-as.data.frame(dBLUPs)  ##!!!!
# Yr20_dBLUPs$Year<-Yr
#   dim(Yr20_dBLUPs)
# Yr20_dBLUPs$plotNo<-expss::vlookup(rownames(Yr20_dBLUPs),dict=dataNHpi[dataNHpi$Year==Yr,],lookup_column = "Crosses",result_column = "plotNo")
# Yr20_dBLUPs$Crosses<-rownames(Yr20_dBLUPs)
# rownames(Yr20_dBLUPs)<-Yr20_dBLUPs$plotNo  
#   Yr20_dBLUPs[Yr20_dBLUPs$Crosses%in%comX,]
#   Yr20_dBLUPs[Yr20_dBLUPs$Crosses=="SL18-LD-13-FG2xSL18-LD-13-MG2",]
# 
# #### Within Yr
# WithinYr_Plot_dBLUPs<-rbind(Yr19_dBLUPs,Yr20_dBLUPs) ##!!!
# save(WithinYr_Plot_dBLUPs,file=paste0(datafdr,"Deregressed_BLUPs_ESplots_plotlevel_WithinYear_AddBD.Rdata"))
#   dim(WithinYr_Plot_dBLUPs)
# 
# WithinYr_Indi_dBLUPs<-rbind(Yr19_dBLUPs,Yr20_dBLUPs)
#   dim(WithinYr_Indi_dBLUPs)
# save(WithinYr_Indi_dBLUPs,file=paste0(datafdr,"Deregressed_BLUPs_ESplots_Indivlevel_WithinYear.Rdata"))

# #### Two Years
# Plot_dBLUPs<-dBLUPs  ##!!!
# save(BLUPs,dBLUPs,file=paste0(datafdr,"Deregressed_BLUPs_ESplots_plotlevel_OverTwoYears_AddBD.Rdata"))
# 
# Indi_dBLUPs<-dBLUPs  ##!!!
# save(BLUPs,dBLUPs,file=paste0(datafdr,"Deregressed_BLUPs_ESplots_Individuallevel_OverTwoYears.Rdata"))
# 
#   dim(Indi_dBLUPs)
#   dim(Plot_dBLUPs)
# sum(rownames(Indi_dBLUPs)%in%rownames(Plot_dBLUPs))
# Both_dBLUPs<-merge(Plot_dBLUPs,Indi_dBLUPs,by.x="row.names",by.y="row.names",all.x=TRUE)  
#   dim(Both_dBLUPs)
#   head(Both_dBLUPs)  
#




## Compare the BLUPs from BothYears vs WithinYear (cor ranged= 0.97-0.98)
both<-merge(Plot_dBLUPs,WithinYr_Plot_dBLUPs,by.x="row.names",by.y="row.names",all.x=TRUE)
for (i in 1:5){
  print(cor(both[,(i+1)],both[,(i+1+5)],use="complete"))
}




#######
#### Summarize the GS TP size--- BothYr
#datafdr<-paste0(WD,"OneTime1920/data/")  ## In terminal
datafdr<-"/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/data/"
load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_Individuals_level_overTwoYears_AddBD.Rdata")) ##!!!

size<-NULL
traits<-colnames(Both_dBLUPs)[-1]
for (t in 1:length(traits)){
  size<-c(size,sum(!is.na(Both_dBLUPs[,traits[t]])|!is.nan(Both_dBLUPs[,traits[t]])))
}
names(size)<-traits
#######
#######

load(paste0(datafdr,"Deregressed_BLUPs_ESplots_plot_Individuals_level_WithinYear_AddBD.Rdata"))


WithinYrSize<-NULL
for (Yr in c(2019,2020)){
  dataf<-droplevels(WithinYr_Both_dBLUPs[WithinYr_Both_dBLUPs$Year.x==Yr,])
  rmcols<-c("Row.names","Year.x","plotNo.x","Crosses.x","Year.y","plotNo.y","Crosses.y")
  traits<-colnames(dataf)[!colnames(dataf)%in%rmcols]
  size<-NULL
  for (t in 1:length(traits)){
    size<-c(size,sum(!is.na(dataf[,traits[t]])|!is.nan(dataf[,traits[t]])))
  }
  names(size)<-traits
  WithinYrSize<-rbind(WithinYrSize,size)
}

