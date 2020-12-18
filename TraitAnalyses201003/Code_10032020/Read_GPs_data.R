
#### Afrer the SNP calling pipeline, 
####!!!!! Run R in terminal, Used "SNPs_Sites.R" script to save csv SNPs data, and to store them in R object for Maggie
###Read in the GPs data
setwd("/Users/maohuang/Desktop/Kelp/2020_2019_Phenotypic_Data/SugarKelpBreeding/TraitAnalyses200820/")
#GPs0<-read.csv("GPsnp_mainGenome.csv",sep=",",header=TRUE,row.names=1)

GPs0<-GPSNP   # GPs0 is the GPSNP/"GPsnp_mainGenome_NA0.8.csv" file

  dim(GPs0)
  GPs0[1:4,1:5]
  GPs0[1631420:nrow(GPs0),183:184]
rownames(GPs0)<-paste0(GPs0$chrom,"_",GPs0$pos)  
map<-GPs0[,1:2]
  head(map)
  tail(map)
GPs<-GPs0[,-c(1:2)]

#### Subset markers for practice purpose
# library(dplyr)
# pick<-sample_n(tbl=map %>% group_by(chrom) ,size=3,replace=TRUE)
#   head(pick)
#   tail(pick)
#   class(pick)  
#   dim(pick)  
# pick<-as.data.frame(pick)
# pick$keep<-paste0(pick$chrom,"_",pick$pos)
# GPs<-GPs0[rownames(GPs0)%in%pick$keep,]
#   dim(GPs)
#   GPs[1:4,1:5]
  

######{{}} Did this in the Terminal  !!!!!!!!!!!!
#### Correct the names of GPs
#write.csv(colnames(GPs),"GPs_sequenced_names.csv")
# GPs_names<-read.csv("GPs_sequenced_names.csv",sep=",",header=T)

# library(stringr)
# names<-str_replace_all(GPs_names[,2],fixed("."),"-")
# names<-str_replace_all(names,"SA-","SA18-")
# names<-str_replace_all(names,"SL-","SL18-")
# names<-str_replace_all(names,"FG-","FG")
# names<-str_replace_all(names,"MG-","MG")

library(stringr)
names<-colnames(GPs)
names<-str_replace_all(names,"SA-","SA18-")
names<-str_replace_all(names,"SL-","SL18-")
names<-str_replace_all(names,"FG-","FG")
names<-str_replace_all(names,"MG-","MG")

colnames(GPs)<-names
  head(names)

######## Impute with A.mat()
# nxm form
tSNP<-t(GPs)
  tSNP[1:5,1:3]
  dim(tSNP)
tSNP[tSNP== 0]=-1   #1st
tSNP[tSNP== 2]=1    #3rd

library(rrBLUP)
K2<-A.mat(tSNP,impute.method="EM",return.imputed=TRUE,shrink=TRUE)    # "A.mat converging:"   # add shrink=TRUE
imputedFromAmat<-K2$imputed
  str(imputedFromAmat)
  imputedFromAmat[1:5,1:6]
AFromAmat<-K2$A    ######!!! 2. GP Amat
hapScoreFromEM <- (imputedFromAmat + 1) / 2  ####
  
GPsA<-AFromAmat


###
M<-imputedFromAmat
GPsA2<-cov(t(M))/mean(diag(cov(t(M))))  #assuming M is (nxp) GenotypesxMarkers allele frequencies matrix (p>>n)

library(ade4)
mantel.rtest(dist(GPsA),dist(GPsA2),nrepet = 999)  #r=0.9732727

save(GPsA,GPsA2,imputedFromAmat,hapScoreFromEM,tSNP,file="GPsAmat_NA0.8.Rdata")  # GPsA2 is estimated by hand


### Look at the common sample #SL18-LD-13-Mg-3 and SL18-LD-13-MG3 are the same sample
cor.test(imputedFromAmat[rownames(imputedFromAmat)=="SL18-LD-13-Mg-3",],imputedFromAmat[rownames(imputedFromAmat)=="SL18-LD-13-MG3",])
#r=0.97
cor.test(tSNP[rownames(tSNP)=="SL18-LD-13-Mg-3",],tSNP[rownames(tSNP)=="SL18-LD-13-MG3",],rm.na=TRUE)
#r=0.97

write.csv(cbind(tSNP[rownames(tSNP)=="SL18-LD-13-Mg-3",],tSNP[rownames(tSNP)=="SL18-LD-13-MG3",]),"Common_sample.csv")  
GPs0common<-cbind(GPs0[,colnames(GPs0)=="SL-LD-13-Mg-3"],GPs0[,colnames(GPs0)=="SL-LD-13-MG-3"])
rownames(GPs0common)<-rownames(GPs0)
write.csv(GPs0common,"Common_sample_from_GPs0.csv")  

#########{{}} Done in Terminal


#library(matrixcal)
source("is.square.matrix.R")
source("is.positive.definite.R")
source("is.symmetric.matrix.R")

load("GPsAmat_NA0.8.Rdata")

# Not Symmetric is due to the decimal points (values not the same for the same sample/off diagonal)

# Further reduced the "SL18-LD-13-Mg-3"
GPsA<-GPsA[!rownames(GPsA)=="SL18-LD-13-Mg-3",!colnames(GPsA)=="SL18-LD-13-Mg-3"]

GPsA<-round(GPsA,digits=5)   ### 2.GPs A
diag(GPsA) <- diag(GPsA) + 1e-5
is.positive.definite(GPsA)  # GPsA2 calculated by hand (above)

fndrsA<-round(mrkRelMat)     ### Added "shrink=TRUE" when estimating the mrkRelMat2
diag(fndrsA) <- diag(fndrsA) + 1e-5
is.positive.definite(fndrsA)   ### 1. fnders A

diag(mrkRelMat2) <- diag(mrkRelMat2) + 1e-5

diag(aMat) <- diag(aMat) + 1e-5
is.positive.definite(aMat)     ### 3. pedigree

save(fndrsA,GPsA,aMat,file="CovList_3_As.Rdata")


### Run these in terminal
load("CovList_3_As.Rdata")
GPsA=newGPsA ############# !!!!!!!! Update once get the new GPsA data

save(fndrsA,GPsA,aMat,file="CovList_3_As.Rdata") ##### !!!! save again

########### Run again in terminal !!!!!!!

load("CovList_3_As.Rdata")
#install.packages("CovCombR")
library(CovCombR)

### 1. NO initial, 3-list
CovList<-NULL
CovList[[1]]<-fndrsA ## fndrsA
CovList[[2]]<-GPsA  ## Further RMed the "SL18-LD-13-Mg-3"
CovList[[3]]<-aMat   
outCovComb1<-CovComb(CovList,nu=1500)   

### 3. amat initial, 3-list
outCovComb3<-CovComb(CovList,nu=1500,Kinit=aMat) 

### 4. amat initial, 3-list, add weight
weights<-c(2,2,1)
outCovComb4<-CovComb(CovList,nu=1500,w=weights,Kinit=aMat) 

save(outCovComb1,outCovComb3,outCovComb4,file="outCovComb_files.Rdata")

# ### 2. aMat initial, 2-list ------- Not correct, Only 239 invididuals
# CovList2<-NULL
# CovList2[[1]]<-fndrsA ## fndrsA
# CovList2[[2]]<-GPsA
# outCovComb2<-CovComb(CovList2,nu=1500,Kinit=aMat) 
####
############ Done in terminal

load("outCovComb_files.Rdata")

#### Then go to the Step2 or Step3 with modeling







####### Not used for now !!!!
#######
#Other Trial in getting the A matrix for the GPs, different ploidy levels
#######
K3<-A.mat(tSNP,impute.method="mean",return.imputed=TRUE)
AFromAmatmean<-K3$A 
hapScoreFromAmatMean <- (K3$imputed + 1) / 2  ####
  hapScoreFromAmatMean[1:5,1:6]
  AFromAmatmean[1:3,1:4]


#### Compare: Impute with mean
imputeToMean <- function(mrkMat){
  imputeVec <- function(mrkVec){
    meanScore <- mean(mrkVec, na.rm=T)
    mrkVec[is.na(mrkVec)] <- meanScore
    return(mrkVec)
  }
  imputedMrkMat <- apply(mrkMat, 2, imputeVec)
  return(imputedMrkMat)
}

# Code as 0,1 and nxm, to use imputeToMean()
tSNP2<-t(GPs)
tSNP2[tSNP2==2]=1
  tSNP2[1:4,1:5]
  t(GPs)[1:4,1:5]
imputedFromMean<-imputeToMean(tSNP2)
  str(imputedFromMean)
  imputedFromMean[1:5,1:6]
  
hapScoreFromMean<-imputedFromMean   ### there are -1??

library(ade4)  
mantel.rtest(dist(imputedFromAmat), dist(imputedFromMean), nrepet = 999) # r=0.979468
mantel.rtest(dist(hapScoreFromEM), dist(hapScoreFromMean), nrepet = 999) # r is the same
# mantel.rtest(dist(hapScoreFromEM), dist(hapScoreFromAmatMean), nrepet = 999) # r=0.979468
# mantel.rtest(dist(hapScoreFromMean), dist(hapScoreFromAmatMean), nrepet = 999) # r=1


### A matrix
 # Calculate GRM as (W %*% W^T) / sum(locusDosageVariance)
  # To calculate p, the function expects dosage coding to be 0, 1, 2 for diploid
  # and 0, 1 for haploid

  calcGenomicRelationshipMatrix <- function(locusMat, ploidy=2){
    if (!any(ploidy == 1:2)) stop("Ploidy must be 1 or 2")
    freq <- colMeans(locusMat) / ploidy
    locusMat <- scale(locusMat, center=T, scale=F)
    return(tcrossprod(locusMat) / sum(ploidy*freq*(1-freq)))
  }
  
FndrCal<-calcGenomicRelationshipMatrix(genoC,ploidy=1) #####???

AFromCal1<-calcGenomicRelationshipMatrix(hapScoreFromEM,ploidy=1)  
AFromCal2<-calcGenomicRelationshipMatrix(hapScoreFromMean,ploidy=1)  

mantel.rtest(dist(AFromAmat*2), dist(AFromCal1), nrepet = 999) # r =0.99 
mantel.rtest(dist(AFromAmat*2), dist(AFromCal2), nrepet = 999) # r= 0.932944

save(GPs,hapScoreFromMean,hapScoreFromEM,AFromAmat,AFromCal1,AFromCal2,file="A.matrix_GPs.Rdata" )


### fndrRows: 79; gpRows:227; spRows: 244
dim(mrkRelMat)
  mrkRelMat[53:56,53:56]
  aMat[53:56,53:56]
dim(AFromAmat)
  AFromAmat[1:4,1:5]
  aMat[79:81,79:81]
  aMat[305:307,305:307]
  
###### All these is to sort he AFromAmat
fndRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) all(vec == 0)))  ### Founders, which() only keeps the TRUE ones
gpRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) is.na(vec[2])))   ### GPs, the col 3 is NA
spRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) all(vec > 0)))    ### Progeny sps, col 2 and 3 are all values
nSp <- length(spRows)
  
  length(fndRows) # 73; 78; 79
  length(gpRows) #99; 156; 227
  length(spRows) #122; 127; 244
  length(fndRows)+length(gpRows)+length(spRows) #294; 361; 550
GP_genotyped_names<-colnames(AFromAmat)  

aMat_GPpart<-aMat[(length(fndRows)+1):(length(fndRows)+length(gpRows)),(length(fndRows)+1):(length(fndRows)+length(gpRows))]
aMat_GPpart_genotyped<-aMat_GPpart[rownames(aMat_GPpart)%in%rownames(AFromAmat),colnames(aMat_GPpart)%in%colnames(AFromAmat)]
#aMat_GPpart_genotyped:105
#AFromAmat: 184

GP_genotyped_ungenotyped_names<-c(intersect(colnames(aMat_GPpart),colnames(AFromAmat)),setdiff(colnames(aMat_GPpart),colnames(AFromAmat)))
AFromAmat_crossed<-AFromAmat[rownames(AFromAmat)%in%intersect(colnames(aMat_GPpart),colnames(AFromAmat)),colnames(AFromAmat)%in%intersect(colnames(aMat_GPpart),colnames(AFromAmat))]
AFromAmat_crossed<-AFromAmat_crossed[intersect(colnames(aMat_GPpart),colnames(AFromAmat)),intersect(colnames(aMat_GPpart),colnames(AFromAmat))]

aMat_GPpart_sort<-aMat_GPpart[GP_genotyped_ungenotyped_names,GP_genotyped_ungenotyped_names]
  dim(aMat_GPpart_sort)
  dim(aMat_GPpart)
  aMat_GPpart_sort[1:5,1:6]
######## How to put aMat_GPpart_sort back to the aMat matrix??? Can aMat be updated so that the genotyped GPs up in the front?  
####### ??????
  
