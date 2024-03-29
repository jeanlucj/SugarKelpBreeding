# Estimate heritabilities with more power (w/o response to E)
# Pedigree tracking for measures of relatedness

# Estimate line and block affect using checks (C-1 and C-2 etc.)
# Estimated breeding value 

# Use called SNPs for marker-based relationship matrix

##########This is JL marker scripting. 
# Output: -1, 0, 1 format matrix with individuals in rows and marker names in columns
# Marker data
# source("convertDArTvcf.R")
# fileName <- "../DArT/Report_DSacc18-3679_SNP_singlerow_2.csv"
# fndrMrkData <- convertDArTvcf(fileName)
# # Filtering of marker data
# # Marker call rate > 0.8
# nNAperMrk <- apply(fndrMrkData, 2, function(v) sum(is.na(v))) / nrow(fndrMrkData)
# fndrMrkData <- fndrMrkData[, nNAperMrk < 0.2]
# nNAperInd <- apply(fndrMrkData, 1, function(v) sum(is.na(v))) / ncol(fndrMrkData)
# # Individual call rate > 0.8
# fndrMrkData <- fndrMrkData[nNAperInd < 0.2,]
# rownames(fndrMrkData) <- paste0(substring(rownames(fndrMrkData), first=1, last=2), "18", substring(rownames(fndrMrkData), first=3))



rm(list=ls())
ls()
library(magrittr)
library(rrBLUP)

# ### Marker Data
# load("Dart2_mNA0.1_nNA0.5_RMhwe0.01_Manual_gl6m_0206.Rdata")
# ### Add 18 to the genotyped founders
# Genotypedname<-paste0(substring(rownames(gl6.m), first=1, last=2), "18", substring(rownames(gl6.m), first=3))
# rownames(gl6.m)<-Genotypedname
# 
#   ls()
#   dim(gl6.m)
#   gl6.m[1:5,1:6]
# 
dataNHpi <-read.csv("dataNHpi_ManuallyAddGrowth_07202020.csv",sep=",",header=TRUE) ##!!!
### RM "checks"
dataNHpi_rmChk<-dataNHpi[!dataNHpi$crossID=="Check",]

FMGPs<-unique(c(as.character(dataNHpi_rmChk$femaPar),as.character(dataNHpi_rmChk$malePar)))
SPs<-strsplit(as.character(FMGPs), split="-", fixed=T)
CrossedSP<-unique(sapply(SPs, function(vec) paste(vec[1:3], collapse="-"))) #93 of them

  head(CrossedSP)
  str(CrossedSP)
  tail(CrossedSP)

# gl6.m2<-gl6.m[rownames(gl6.m)%in%CrossedSP,]
#   dim(gl6.m2)

load("FarmCPU_GAPIT.Rdata")
#write.csv(geno2$taxa,"Fndr_Crossed_genotyped_Dart1.csv")
  ls()
geno2<-geno[,-1]
rownames(geno2)<-geno$taxa
  geno2[1:5,1:6]
  dim(geno2)
rownames(geno2) <- paste0(substring(rownames(geno2), first=1, last=2), "18", substring(rownames(geno2), first=3))
geno2[geno2==0]=-1
geno2[geno2==1]=0
geno2[geno2==2]=1

geno2<-geno2[rownames(geno2)%in%CrossedSP,]  # 93 fnders made crosses, only 56 were genotyped
  

#write.csv(rownames(geno2),"Fndrs_genotyped_Samples.csv")
fndrMrkData<-geno2
mrkRelMat <- A.mat(fndrMrkData, impute.method="EM", return.imputed=T)
fndrMrkDataImp <- mrkRelMat$imputed
mrkRelMat <- mrkRelMat$A
  dim(mrkRelMat)  #56 x 56
#write.csv(rownames(mrkRelMat),"Fndrs_used_in_Amat.csv")
  
########################### A. Plot
dataNHpi <-read.csv("dataNHpi_ManuallyAddGrowth_07202020.csv",sep=",",header=TRUE) ##!!!
dataNHpi <- dataNHpi[order(dataNHpi$plotNo),]  #### Plots in alphabetic order

  dim(dataNHpi)
  str(dataNHpi)
dataNHpi<-dataNHpi[dataNHpi$Region=="GOM",] ## !!! RM SNE, 530 rows
    
dataNHpi<-dataNHpi[!dataNHpi$crossID=="Buffer",]  ## !!! RMed Buffer lines

#dataNHpi<-dataNHpi[!dataNHpi$PhotoScore==0,] ## !!! RM PhotoScore=0,  447 rows
dataNHpi<-dataNHpi[dataNHpi$PhotoScore>1,] ## !!! RM PhotoScore<1

  dim(dataNHpi)
  colnames(dataNHpi)

dataNHpi19<-dataNHpi[dataNHpi$Year==2019,] ## 2019 data   !!!!!!
dataNHpi19_C<-dataNHpi19
dataNHpi19_C<-dataNHpi19_C[order(dataNHpi19_C$plotNo),]  ## Order plotNo alphabetically

dataNHpi20<-dataNHpi[dataNHpi$Year==2020,]   #
dataNHpi20_C<-dataNHpi20
dataNHpi20_C<-dataNHpi20_C[order(dataNHpi20_C$plotNo),]  ## Order plotNo alphabetically

dataNHpiBoth_C<-dataNHpi  ## Order plotNo alphabetically
dataNHpiBoth_C<-dataNHpiBoth_C[order(dataNHpiBoth_C$plotNo),] ## Order plotNo alphabetically


save(dataNHpi19_C,dataNHpi20_C,dataNHpiBoth_C,file="dataNHpi_withChk_3_sets_PhotoScore23.rdata")

## RM checks to make the pedigree-biphasic matrix and CC-matrix
dataNHpi19<-dataNHpi19[!dataNHpi19$crossID=="Check",]  
dataNHpi20<-dataNHpi20[!dataNHpi20$crossID=="Check",]  
dataNHpiBoth<-dataNHpi[!dataNHpi$crossID=="Check",]

  dim(dataNHpi19) #122
  dim(dataNHpi20)  #128
  dim(dataNHpiBoth)  #250

######################################################################  
### Make pedigree relationship matrix  !!!!!!!!!! RUN 3 TIMES !!!!!! 

#1. Create the pedigree
source("makeBiphasicPed.R")

dataNHpi<-dataNHpi19  ######### !!!!
yr<-"2019"


dataNHpi<-dataNHpi20  ######### !!!!
yr<-"2020"


dataNHpi<-dataNHpiBoth ########## !!!!!
yr<-"Both"


dataNHpi$Crosses<-as.factor(as.character(dataNHpi$Crosses))
dataNHpi$femaPar<-as.character(dataNHpi$femaPar)
dataNHpi$malePar<-as.character(dataNHpi$malePar)

dataNHpi<-dataNHpi[order(dataNHpi$plotNo),]  ### Order plots alphabetically, so that order match that in making the Z matrix

kelpNameColumns <- dataNHpi[, c("femaPar", "malePar")]
kelpNameColumns <- kelpNameColumns[kelpNameColumns[,1] != "",]   #### ? Is this removing the Checks?

### Need to remove the femaPar x malePar that's the same cross
kelpNameColumns$RM<-paste0(kelpNameColumns[,1],"x",kelpNameColumns[,2])
kelpNameColumns<-kelpNameColumns[!duplicated(kelpNameColumns$RM),]
kelpNameColumns<-kelpNameColumns[,-3]

biphasicPedNH <- makeBiphasicPed(kelpNameColumns, rownames(mrkRelMat)) #### The rownames(mrkRelMat) is to put the fndrs in the front in the pedigree matrix
  head(kelpNameColumns)
  head(rownames(mrkRelMat))

write.csv(biphasicPedNH,paste0("biphasicPedNH_",yr,"_PhotoScore23.csv") ) 
      

# 2.Calculate the relationship matrix
source("calcCCmatrixBiphasic.R")
biphasicCCmat <- calcCCmatrixBiphasic(biphasicPedNH)  #### !!!!!!! Update into calcCCmatrixHaploid() ?????

rownames(biphasicCCmat) <- colnames(biphasicCCmat) <- rownames(biphasicPedNH)


# 3. Combine the marker- with the pedigree- relationship matrix to make H matrix #### !!!!! How to do that?????
aMat <- 2 * biphasicCCmat
  dim(aMat)
fndRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) all(vec == 0)))  ### Founders, which() only keeps the TRUE ones
gpRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) is.na(vec[2])))   ### GPs, the col 3 is NA
spRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) all(vec > 0)))    ### Progeny sps, col 2 and 3 are all values
nSp <- length(spRows)

  length(fndRows) # 73; 78; 79
  length(gpRows) #99; 156; 227
  length(spRows) #122; 127; 244
  length(fndRows)+length(gpRows)+length(spRows) #294; 361; 550
hMat <- calcHmatrix(mrkRelMat, aMat, aMatFounders=rownames(mrkRelMat))

save(mrkRelMat,aMat,hMat,biphasicPedNH,biphasicCCmat,fndrMrkData, file=paste0("hMat_PedNH_CCmat_fndrMrkData_",yr,"_PhotoScore23.rdata"))  ###### !!!!!!!



############# B. Individual
dataNHim <- read.csv("Indi_2019_2020_Compile_07202020_Edit.csv",sep=",",header=TRUE)
  dim(dataNHim)
  str(dataNHim)
  tail(dataNHim)

dataNHim<-dataNHim[dataNHim$Region=="GOM",] ## !!! RM SNE, 530 rows
  dim(dataNHim)  #5046

dataNHim<-dataNHim[!dataNHim$crossID=="Buffer",]
  dim(dataNHim) #4996

dataNHim19_C<-dataNHim[dataNHim$Year==2019,] 
dataNHim19_C<-dataNHim19_C[order(dataNHim19_C$plotNo),]  ### Plots in Alphabetic order; Here the plot did not RM plots with photoScore < 2,3
  dim(dataNHim19_C)  # 2800 x 12
  tail(dataNHim19_C)

dataNHim20_C<-dataNHim[dataNHim$Year==2020,]
dataNHim20_C<-dataNHim20_C[order(dataNHim20_C$plotNo),]
  dim(dataNHim20_C)  # 2550 x 12
  tail(dataNHim20_C)
  
dataNHimboth_C<-dataNHim
dataNHimboth_C<-dataNHimboth_C[order(dataNHimboth_C$plotNo),]
  tail(dataNHimboth_C)
###dataNHim$plotNo <- gsub("S", "", dataNHim$plotNo, fixed=TRUE) ### The 2020 plotNo is "SXXX"; 2019 is just number
###dataNHim$plotNo[dataNHim$plotNo == ""] <- "C2-I"



save(dataNHim19_C,dataNHim20_C,dataNHimboth_C,file="dataNHim_withChk_3_sets_PhotoScore0123.rdata")
#### !!!! This individual dataset still HAS the photoscore < 2,3 plots




