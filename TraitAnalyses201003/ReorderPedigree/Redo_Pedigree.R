# Re do the pedigree matrix
save(GPsequenced,kelpNameColumns,SPGPs,file="What_has_been_Used_from_Step1_Create_Matrix.Rdata")

##### Re-order the pedigree matrix
setwd("/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/TraitAnalyses201003/Making_haploid_CovComb/InputData")
biphasicPedNH<-read.csv("biphasiPedNH_addmoreGP_866.csv",sep=",",header=TRUE,row.names=1)


###1 GP genotyped,crossed and had phenotypic data in dataNHpi
###2 GP not genotyped?,crossed and had phenotypic data in dataNHpi

###3 GP genotyped,not crossed and no phenotypic data in dataNHpi
###4 GP from 2019S_plot level

#### "GPsequenced" has the list of GPs being sequenced in 3 plates
#### "kelpNameColumns" has the list of unique femaPar x malePar, GPs made cross and had phenotypic data
#### SPGPs has the list of SP_GPs to be included

str(GPsequenced)  # 270
str(kelpNameColumns)  # 244
str(unique(GPsequenced)) #270
str(unique(kelpNameColumns$femaPar))  #106
str(unique(kelpNameColumns$malePar))  #121

##1 
FG_seq<-unique(kelpNameColumns$femaPar[kelpNameColumns$femaPar%in%GPsequenced]) #68
str(FG_seq)
MG_seq<-unique(kelpNameColumns$malePar[kelpNameColumns$malePar%in%GPsequenced]) #91
str(MG_seq)

GP_seq<-c(FG_seq,MG_seq)   # a Unique list for GPs in 1
##2
FG_Noseq<-unique(kelpNameColumns$femaPar[!kelpNameColumns$femaPar%in%GPsequenced]) #38
str(FG_Noseq)
MG_Noseq<-unique(kelpNameColumns$malePar[!kelpNameColumns$malePar%in%GPsequenced]) #30
str(MG_Noseq)

GP_Noseq<-c(FG_Noseq,MG_Noseq)  # a Unique list for GPs in 2

##3  
GP_noPheno<-GPsequenced[!GPsequenced%in%c(kelpNameColumns$femaPar,kelpNameColumns$malePar)] #111
str(GP_noPheno)             # a Unique list for GPs in 3

##4 
S_GP<-SPGPs$GametophyteID


# I. find all 1,2,3,4 their fnders, order, unique the list
# II. GPs themselves, 1,2 and 3
# III.The pair of SP plots in kelpNamesColumns femaParx malePar
# IV. The rest of the SPGPs$SPcrosses (parents for 4); adding that plot for SPGPs$SPCrosses that is not in the kelpNamesColumns femaPar x malePar
# V. S_GPs for 4 themselves

find_fnder_for_GP<-function(GPstring){
  fndr<-strsplit(GPstring,split="-", fixed=T) 
  fnder_unique<-unique(sapply(fndr,function(vec) paste(vec[1:3], collapse="-")))
  return(fnder_unique)
}

#1
GP_seq_fndr<-find_fnder_for_GP(GP_seq) #64
#2
GP_Noseq_fndr<-find_fnder_for_GP(GP_Noseq) # 44
str(GP_seq_fndr)  # 1
str(GP_Noseq_fndr) # 2
str(unique(c(GP_seq_fndr,GP_Noseq_fndr)))  # 70 
#3
GP_noPheno_fndr<-find_fnder_for_GP(GP_noPheno)  #61 
str(GP_noPheno_fndr) # 3
str(unique(c(unique(c(GP_seq_fndr,GP_Noseq_fndr)),GP_noPheno_fndr)))  #93

library(stringr)
find_fnder_for_SGP<-function(S_GPCross){
  F_M_GPs<-str_split_fixed(string=as.character(S_GPCross), "x", 2)
  S_FGP<-F_M_GPs[,1]
  S_MGP<-F_M_GPs[,2]
  S_FGP_fndr<-find_fnder_for_GP(S_FGP)  #11
  S_MGP_fndr<-find_fnder_for_GP(S_MGP)  #9
  
  return(list(S_FGP=S_FGP,S_MGP=S_MGP,S_FGP_fndr=S_FGP_fndr,S_MGP_fndr=S_MGP_fndr))
}


S_GP_its_fnders<-c(find_fnder_for_SGP(SPGPs$SPCross)$S_FGP_fndr,find_fnder_for_SGP(SPGPs$SPCross)$S_MGP_fndr)
str(S_GP_its_fnders)
str(unique(c(unique(c(unique(c(GP_seq_fndr,GP_Noseq_fndr)),GP_noPheno_fndr)),S_GP_its_fnders))) # 93 

S_GP_its_parentalGP<-unique(c(find_fnder_for_SGP(SPGPs$SPCross)$S_FGP,find_fnder_for_SGP(SPGPs$SPCross)$S_MGP))
str(S_GP_its_parentalGP)

#For 4, all S_GP_its_fnders were already included
#For 1, to row64, For 2 row65 to row70, For 3 row71 to row93)
#### List of fnders for 1,2,3,4

All_fnders<-c(GP_seq_fndr,GP_Noseq_fndr,GP_noPheno_fndr,S_GP_its_fnders)  # 189
str(All_fnders)
All_fnders_unique<-unique(All_fnders)   #93
str(All_fnders_unique)

## I. Ped for fnders
nFounders<-length(All_fnders_unique) 
fndRow <- 1:nFounders
names(fndRow) <- All_fnders_unique
fndPed <- cbind(fndRow, 0, 0)   


All_GPs123<-c(GP_seq,GP_Noseq,GP_noPheno)
str(All_GPs123)     # 338, unique 
length(GP_seq)    # row94 to row252 (94+159-1)
length(GP_Noseq)  # row253 to row320 (253+68-1)
length(GP_noPheno) #row 321 to row431 (321+111-1)
sum(S_GP_its_parentalGP %in%(All_GPs123))==length(S_GP_its_parentalGP)   #All the GPs used for 4 are in the All_GPs123 list

## II. Ped for GPs #### DO I NEED TO ORDER THEM BY FEMALE firt and then BY MALE ???
GPsRow <- nFounders + 1:length(All_GPs123)
names(GPsRow) <- All_GPs123
founderPar<-sapply(strsplit(All_GPs123, split="-", fixed=T), function(vec) paste(vec[1:3], collapse="-"))
GPsPed <- cbind(GPsRow, fndRow[founderPar], NA)

## III. Ped for the SP plots in kelpNamesColumns femaParx malePar          
sporProg1 <- 1:nrow(kelpNameColumns) + nrow(fndPed) + nrow(GPsPed)
names(sporProg1)<-paste0(kelpNameColumns[,1],"x",kelpNameColumns[,2])
progPed1 <- cbind(sporProg1, GPsRow[kelpNameColumns[,1]], GPsRow[kelpNameColumns[,2]]) # row 432 to row 675 (431+244)           

## IV. Ped for the SP_GP's plots that were not in 3        
# The SPplots not listed in the kelpNameColumns combination (because their plot photo score may not pass 1)
sproProg2Name<-unique(droplevels(SPGPs$SPCross[!SPGPs$SPCross%in%names(sporProg1)]))        
sproProg2<-1:length(sproProg2Name)+nrow(fndPed)+nrow(GPsPed)+nrow(progPed1)

names(sproProg2)<-sproProg2Name
progPed2<-cbind(sproProg2,GPsRow[str_split_fixed(string=as.character(names(sproProg2)), "x", 2)[,1]], GPsRow[str_split_fixed(string=as.character(names(sproProg2)), "x", 2)[,2]])

progPed<-rbind(progPed1,progPed2)

## V. Ped for the SP_GPs themselves
spRows<-progPed[,1]
nrow(SPGPs)==length(unique(SPGPs$GametophyteID))
SP_GP_row<-1:nrow(SPGPs)+nrow(fndPed)+nrow(GPsPed)+nrow(progPed)
names(SP_GP_row)<-SPGPs$GametophyteID
SP_GP_Ped<-cbind(SP_GP_row,spRows[as.character(SPGPs$SPCross)],NA)
head(SP_GP_Ped)

Ped_in_Order<-rbind(fndPed,GPsPed,progPed,SP_GP_Ped)  
dim(Ped_in_Order)
write.csv(Ped_in_Order,"Ped_in_Order_754_Individuals.csv")


