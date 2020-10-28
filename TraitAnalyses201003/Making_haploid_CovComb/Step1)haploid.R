
###### This is to work the  the haploid level of matrices
#### 1. Founders marker data, make haploid level, one sample 2 rows: 0->0,0 1->0,1 2->1,1
#### Cal relationship matrix for the fnders
#### 2. Cal hapccMatrix for the pedigree
#### 3. GPs relationship matrix, 0->0,2->1
#### 4. Use CovComb() to combine this Founder relationship amat, and the GP sequenced amat and the PedNH ccmatrix
#### 5. condense to diploid level
#######

#### function to make diploid marker data to be haploid
#!!!!!! This is for haploid
######### Diploid split into 2 rows (onerow at a time) format
ConvertOneRow<-function(fndrMrkDataImp){
onerow<-fndrMrkDataImp
onerow[onerow< -1]=-1
onerow[onerow> -1 & onerow< -0.5]=-1
onerow[onerow> -0.5 & onerow< 0]=0
onerow[onerow> 0 & onerow< 0.5]=0
onerow[onerow> 0.5 & onerow< 1]=1
onerow[onerow> 1]=1
### Two row format for haploid
geno3<-onerow
geno3[geno3==-1]=0

geno4<-onerow
geno4[geno4==1]=1
geno4[geno4==0]=1
geno4[geno4==-1]=0
  geno[1:4,1:5]
  fndrMrkDataImp[1:4,1:5]
  onerow[1:4,1:5]
  geno3[1:4,1:5]
  geno4[1:4,1:5]
geno4<-as.matrix(geno4)
geno3<-as.matrix(geno3)
genoC<-matrix(nrow=nrow(geno4)*2,ncol=ncol(geno4))
colnames(genoC)<-colnames(geno4)

for (i in 1:nrow(geno4)){
  genoC[((2*i)-1),]=as.vector(geno3[i,])
  genoC[(2*i),]=geno4[i,]
  genoC<-as.data.frame(genoC)
  rownames(genoC)[(2*i)-1]<-rownames(geno3)[i]
  rownames(genoC)[(2*i)]<-paste0(rownames(geno4),"_2")[i]
  }

return(list(OneRowG=genoC))
}
#####

#### 1. fnder marker data
load("FarmCPU_GAPIT.Rdata")
ls()
geno2<-geno[,-1]
rownames(geno2)<-geno$taxa
geno2[1:5,1:6]
dim(geno2)
rownames(geno2) <- paste0(substring(rownames(geno2), first=1, last=2), "18", substring(rownames(geno2), first=3))
geno2[geno2==0]=-1
geno2[geno2==1]=0
geno2[geno2==2]=1

#######Did not do this step, 
#so some fnders were not in the pedigree file but were included to estimate aMat, they were removed out later
#though this should not be an issue??
#######
#geno2<-geno2[rownames(geno2)%in%CrossedSP,]  # 93 fnders made crosses, only 56 were genotyped
#write.csv(rownames(geno2),"Fndrs_genotyped_Samples.csv")

fndrMrkData<-geno2

## Impute fndr Mrk at diploid level
mrkRelMat <- A.mat(fndrMrkData, impute.method="EM", return.imputed=T,shrink=TRUE) ## Add shrink per Deniz
fndrMrkDataImp <- mrkRelMat$imputed

fndOneRow<-ConvertOneRow(fndrMrkDataImp)$OneRowG
  dim(fndOneRow)
  fndOneRow[1:4,1:5]  
  dim(fndrMrkDataImp)

# Calculate GRM as (W %*% W^T) / sum(locusDosageVariance)
# To calculate p, the function expects dosage coding to be 0, 1, 2 for diploid
# and 0, 1 for haploid
  calcGenomicRelationshipMatrix <- function(locusMat, ploidy=2){
    if (!any(ploidy == 1:2)) stop("Ploidy must be 1 or 2")
    freq <- colMeans(locusMat) / ploidy
    locusMat <- scale(locusMat, center=T, scale=F)
    return(tcrossprod(locusMat) / sum(ploidy*freq*(1-freq)))
  }
  
fndrA<-calcGenomicRelationshipMatrix(fndOneRow,ploidy=1)  ###!!! 1. fnder marker Rel matrix, haploid
#fndrAcondese<-condenseHapCCmat(fndrA)

##Another way of estimating the relationship Matrix. Why results not the same??
#M<-fndOneRow  ##### RelationshipMatrix for GP another way
#fndrA2<-cov(t(M))/mean(diag(cov(t(M))))  #assuming M is (nxp) GenotypesxMarkers allele frequencies matrix (p>>n)

  
## GPs change rows format    

#### 2. ####  
setwd("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/Phenotypic_Analysis/TraitAnalyses200820_Updated_AfterCrossList/withSGP")
biphasicPedNH<-read.csv("biphasiPedNH_addmoreGP_866.csv",sep=",",header=T,row.names=1) 
  dim(biphasicPedNH)
  head(biphasicPedNH)
  tail(biphasicPedNH)
# #########
  
source("calcCCmatrixBiphasic.R")  
  
biphasicPedNH<-as.matrix(biphasicPedNH)  
#namesNum<-grep("UCONN",rownames(biphasicPedNH))
#namesLS<-rownames(biphasicPedNH)[namesNum]
#library(stringr)
#rownames(biphasicPedNH)[namesNum]<-apply(str_split_fixed(namesLS,"-",4), 1, function(row) paste(row[c(1,3,4)], collapse="-"))

HaploidCCcal <- calcCCmatrixHaploid(biphasicPedNH)

biphasichapccMat<-HaploidCCcal$hccMat   ###!!! 2. CCmatrix/aMat at haploid level for pedigree
hapOnePointer<-HaploidCCcal$hapOnePointer  ###!!! For condensing
dipTwoPointers<-HaploidCCcal$dipTwoPointers  ###!!! For condensing

### This is to rename the rownames out of the biphasichapccMat
pedinames<-NULL
for (i in 1:nrow(biphasicPedNH)){
  if(!is.na(biphasicPedNH[i,2]) & !is.na(biphasicPedNH[i,3])){
    # this is diploid, both parents have numbers
    pedinames<-c(pedinames,rownames(biphasicPedNH)[i],paste0(rownames(biphasicPedNH)[i],"_2"))

  }else{
    # this is GP, one parent is missing
    pedinames<-c(pedinames,rownames(biphasicPedNH)[i])
  }
}

rownames(biphasichapccMat)<-colnames(biphasichapccMat)<-pedinames

biphasicmCCmat<-HaploidCCcal$mccMat ### This is the condensed hapCC matrix
# This is identical to the biphasicmCCmat as it use the same function
# biphasichapccMat_condense<-condenseHapCCmat(biphasichapccMat,hapOnePointer,dipTwoPointers)

load("hMat_PedNH_CCmat_fndrMrkData_Both_PhotoScore23_WithSGP_866.rdata")
biphasciCCmat<-aMat/2 # This is the CCmatrix estimated at diploid level
library(cultevo)
mantel.test(dist(as.matrix(biphasicmCCmat)),dist(as.matrix(aMat)),trials=99) # r=1
mantel.test(dist(as.matrix(biphasicCCmat)),dist(as.matrix(biphasichapccMat_condense)),trials=99) # r=1

biphasichapccMat[140:142,140:142]

#### 3. GPs data
load("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/GenotypicData_for_SugarKelpBreeding/GPs_mainGenome_NA0.8_P1P2P3.Rdata")
# GPSNP is raw SNPs

load("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/GenotypicData_Files_for_Phenotypic_Analysis/GPsAmat_NA0.8_P1P2P3_09282020.Rdata")
  dim(imputedFromAmat)
  imputedFromAmat[1:4,1:5]

#### Need to redo the formating of the raw marker data again  
  GPs0<-GPSNP   # 
  rownames(GPs0)<-paste0(GPs0$chrom,"_",GPs0$pos)  
  map<-GPs0[,1:2]
  GPs<-GPs0[,-c(1:2)]
  
  ######{{}} Did this in the Terminal to estimate the RelationshipMatrix2 
  #### Correct the names of GPs
  
  library(stringr)
  names<-colnames(GPs)
  names<-str_replace_all(names,"SA-","SA18-")
  names<-str_replace_all(names,"SL-","SL18-")
  names<-str_replace_all(names,"FG-","FG")
  names<-str_replace_all(names,"MG-","MG")
  
  colnames(GPs)<-names
  head(names)
  
  ## RM checks, the mistake one (RM SL-PI-1-FG-1,SL-CT1-FG-3)
  IndiRM<-c("SL18-PI-1-FG1","3","SL18-CT1-FG3","SL18-CT1-MG2","SL18-CT1-MG3","SL18-OI-15-Female","SL18-ME-1-FG1","SL18-ME-1-MG1","SL18-ME-1-MG2")
  GPs1<-GPs[,!colnames(GPs)%in%IndiRM]
  GPs<-GPs1
  
  #E. Impute GP marker matrix with A.mat()
  tSNP<-t(GPs)
    dim(tSNP)
  
  #RMlist<-c("SL18-PI-1-FG1","SL18-CT1-FG3")
  #GPsMrk<-imputedFromAmat[-which(rownames(imputedFromAmat)%in%RMlist),]
  #dim(GPsMrk)

  tSNP[tSNP==2]=1  # Change to haploid level dosage
    
#mrkMat needs to be m xn dimension
  
  imputeToMean <- function(mrkMat){
    imputeVec <- function(mrkVec){
      meanScore <- mean(mrkVec, na.rm=T)
      mrkVec[is.na(mrkVec)] <- meanScore
      return(mrkVec)
    }
    imputedMrkMat <- apply(mrkMat, 2, imputeVec)
    return(imputedMrkMat)
  }
  
GPsMrkImp<-imputeToMean(tSNP)  
  dim(GPsMrkImp)
  GPsMrkImp[1:4,1:5]
  tSNP[1:4,1:5]
  
onerow<-GPsMrkImp 
  onerow[onerow> 0 & onerow< 0.5]=0
  onerow[onerow>= 0.5 & onerow< 1]=1
  onerow[onerow> 1]=1  

GPsA<-calcGenomicRelationshipMatrix(GPsMrkImp,ploidy=1)
GPsA2<-calcGenomicRelationshipMatrix(onerow,ploidy=1) ### This function assumes data to be 0 and 1 for haploid

mantel.test(dist(GPsA),dist(GPsA2),trials=99) # 0.969
  
M<-GPsMrkImp
GPsA3<-cov(t(M))/mean(diag(cov(t(M))))
mantel.test(dist(GPsA),dist(GPsA3),trials=99) # 0.984

M2<-onerow
GPsA4<-cov(t(M2))/mean(diag(cov(t(M))))
mantel.test(dist(GPsA2),dist(GPsA4),trials=99) # 0.9

### Choose the onerow, calcGenomicRelationshipMatrix(), so GPsA2

### This is the GPsA2

### 4. ###
library(CovCombR)
source("is.square.matrix.R")
source("is.positive.definite.R")
source("is.symmetric.matrix.R")

diag(fndrA) <- diag(fndrA) + 1e-5
is.positive.definite(fndrA)
diag(GPsA) <- diag(GPsA) + 1e-5
is.positive.definite(GPsA)
diag(biphasichapccMat) <- diag(biphasichapccMat) + 1e-5
is.positive.definite(biphasichapccMat)

### RM the extra fndr individuals
RMextraFndr<-which(!rownames(fndrA)%in%rownames(biphasichapccMat))
fndrA2<-fndrA[-RMextraFndr,-RMextraFndr]
  dim(fndrA2)
  fndrA2[29:30,29:30]
### Is GPsA2 individuals all in the biphasichapccMat?
which(!rownames(GPsA2)%in%rownames(biphasichapccMat))  
is.positive.definite(fndrA2)  
save(fndrA2,GPsA2,biphasichapccMat,file="Newly_saved_3_As_for_CovComb.Rdata")
save(hapOnePointer,dipTwoPointers,file="Pointers_from_biphasichapccMat.Rdata")
##{{Done in terminal}}
CovList<-NULL
CovList[[1]]<-fndrA2 ## fndrsA
  
CovList[[2]]<-GPsA2  
CovList[[3]]<- biphasichapccMat  

weights<-c(2,2,1)
outCovComb4<-CovComb(Klist=CovList,nu=2000,w=weights,Kinit=biphasichapccMat) 

source("calcCCmatrixBiphasic.R")
outCovComb4_conden<-condenseHapCCmat(outCovComb4,hapOnePointer,dipTwoPointers)

save(outCovComb4,outCovComb4_conden,file="outCovComb4_and_Conden.Rdata")
###{{}}



# run source code
Kinit=biphasichapccMat
Klist<-CovList
Kinvlist=NULL
Kinvlist<-lapply(Klist,solve)
  dim(Kinvlist)

linenames <- c()    
for (i in seq_len(length(Kinvlist))){
    linenames <- union(linenames, rownames(Kinvlist[[i]]))
  }

Kinit <- Kinit[match(linenames, colnames(Kinit)), match(linenames, 
                                                        rownames(Kinit))]/nu
  

#save(fndrA,GPsA2,biphasichapccMat,file="All_3_RelationshipMatrices.Rdata")
