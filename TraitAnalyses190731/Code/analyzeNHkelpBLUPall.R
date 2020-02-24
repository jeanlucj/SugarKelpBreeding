# NOTE: Set working directory to that containing this code [called "Code" here]
# There should be a directory called "20190719" at the same level as "Code"
#  which has the phenotypic data
# There should be a directory called "DArT" which has the genetic data
library(magrittr)
library(rrBLUP)
# spRows has 199 elements. dataNHpi has 218 rows.  The last 19 rows are checks
dataNHpi <- read.csv("../20190719/NH Phenotyping Data - Plot Info.csv", header=T, as.is=T)
dataNHpi <- dataNHpi[,c(1, 3:6, 8:10, 14:15, 28:29, 42:44, 47:51)]
colnames(dataNHpi) <- c("plotNo", "femaPar", "femaParLoc", "malePar", "maleParLoc", "line", "block", "date", "wetWgtPlot", "lengthPlot", "numBlades5cm", "densityBlades", "percDryWgt", "wetWgtPerM", "dryWgtPerM", "moistPerc", "protPerc", "fatPerc", "fiberPerc", "ashPerc")

dataNHim <- read.csv("../20190719/NH Phenotyping Data - Individual Measurements.csv", header=T, as.is=T)
dataNHim <- dataNHim[, c(1, 3:12)]
colnames(dataNHim) <- c("plotNo", "subPlot", "bladeLength", "bladeMaxWidth", "bladeWidth10cm", "bladeThickness", "stipeLength", "stipeDiameter", "stipeHollow", "bullations", "fouling")
dataNHim$plotNo <- gsub("S-", "", dataNHim$plotNo, fixed=TRUE)
dataNHim$plotNo[dataNHim$plotNo == ""] <- "C2-I"

source("convertDArTvcf.R")
fileName <- "../DArT/Report_DSacc18-3679_SNP_singlerow_2.csv"
fndrMrkData <- convertDArTvcf(fileName)
# Filtering of marker data
nNAperMrk <- apply(fndrMrkData, 2, function(v) sum(is.na(v))) / nrow(fndrMrkData)
fndrMrkData <- fndrMrkData[, nNAperMrk < 0.2]
nNAperInd <- apply(fndrMrkData, 1, function(v) sum(is.na(v))) / ncol(fndrMrkData)
fndrMrkData <- fndrMrkData[nNAperInd < 0.2,]
rownames(fndrMrkData) <- paste0(substring(rownames(fndrMrkData), first=1, last=2), "18", substring(rownames(fndrMrkData), first=3))
mrkRelMat <- A.mat(fndrMrkData, impute.method="EM", return.imputed=T, n.core=4)
fndrMrkDataImp <- mrkRelMat$imputed
mrkRelMat <- mrkRelMat$A

source("makeBiphasicPed.R")
kelpNameColumns <- dataNHpi[, c("femaPar", "malePar")]
kelpNameColumns <- kelpNameColumns[kelpNameColumns[,1] != "",]
biphasicPedNH <- makeBiphasicPed(kelpNameColumns, rownames(mrkRelMat))

source("calcCCmatrixBiphasic.R")
biphasicCCmat <- calcCCmatrixBiphasic(biphasicPedNH)
rownames(biphasicCCmat) <- colnames(biphasicCCmat) <- rownames(biphasicPedNH)
aMat <- 2 * biphasicCCmat
fndRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) all(vec == 0)))
gpRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) is.na(vec[2])))
spRows <- which(apply(biphasicPedNH[,2:3], 1, function(vec) all(vec > 0)))
nSp <- length(spRows)

hMat <- calcHmatrix(mrkRelMat, aMat, aMatFounders=rownames(mrkRelMat))

###################################
## Data analysis really starts here
# Plot information analysis
dataNHpi$popChk <- ifelse(substr(dataNHpi$plotNo, 1, 1) == "C", substr(dataNHpi$plotNo, 1, 2), "ES")
dataNHpi$withinLoc <- c(ifelse(dataNHpi$femaParLoc[1:nSp] == dataNHpi$maleParLoc[1:nSp], 1, 0), rep(0, 19))

for (col in c("wetWgtPlot", "lengthPlot", "numBlades5cm", "densityBlades", "percDryWgt", "wetWgtPerM", "dryWgtPerM")) dataNHpi[,col] <- as.numeric(dataNHpi[,col])
dataNHpi[is.na(dataNHpi$percDryWgt), c("wetWgtPerM", "dryWgtPerM")] <- NA
keepRows <- !is.na(dataNHpi$percDryWgt)

isSelf <- integer(nrow(dataNHpi))
isSelf[diag(aMat[spRows, spRows]) > 1] <- 1
dataNHpi$isSelf <- isSelf

for (col in c("plotNo", "femaPar", "femaParLoc", "malePar", "maleParLoc", "line", "block", "date", "popChk", "withinLoc", "isSelf")) dataNHpi[,col] <- as.factor(dataNHpi[,col])

msX <- model.matrix( ~ line + block + date + popChk, data=dataNHpi)
# Construction of Z assumes that all checks are at the bottom of the dataset
msZ <- rbind(cbind(matrix(0, nSp, nrow(aMat) - nSp), diag(nSp)), matrix(0, nrow(dataNHpi) - nSp, nrow(aMat)))

msOutWWP1 <- mixed.solve(y=log(dataNHpi$wetWgtPlot+1), Z=msZ, K=aMat, X=msX)
msOutDB1 <- mixed.solve(y=dataNHpi$densityBlades, Z=msZ, K=aMat, X=msX)
msOutDWPM1 <- mixed.solve(y=log(dataNHpi$dryWgtPerM+1), Z=msZ, K=aMat, X=msX)
msOutPDW1 <- mixed.solve(y=dataNHpi$percDryWgt, Z=msZ, K=aMat, X=msX)
heritability <- function(msOut){
  return(msOut$Vu / (msOut$Vu + msOut$Ve))
}
h2naive <- c(heritability(msOutWWP1), heritability(msOutDB1), heritability(msOutDWPM1), heritability(msOutPDW1))
names(h2naive) <- c("wetWgtPlot", "densityBlades", "dryWgtPerM", "percDryWgt")

# Add in two components:
# 1. Density Blades covariate
# 2. Replace fixed effect of experimental versus check with the experimental location of origin
# Both components will take out some genetic variation.  Should probably figure out how much each...
femaLocInc <- model.matrix(~ -1 + femaParLoc, data=dataNHpi)
maleLocInc <- model.matrix(~ -1 + maleParLoc, data=dataNHpi)
locInc <- (femaLocInc[,-1] + maleLocInc[,-1]) / 2

msXdb <- model.matrix( ~ densityBlades + line + block + date + popChk, data=dataNHpi)
# locInc sums to 1 for all experimental sporophytes but 0 for the checks
# So perfectly colinear with the popChkES incidence column (17)
# msXdb <- cbind(msX[, -17], locInc) # 20190720 locInc causing a problem.  Not using it.

msOutWWP <- mixed.solve(y=log(dataNHpi$wetWgtPlot+1), Z=msZ, K=aMat, X=msXdb, SE=T)
msOutDWPM <- mixed.solve(y=log(dataNHpi$dryWgtPerM+1), Z=msZ, K=aMat, X=msXdb, SE=T)
msOutPDW <- mixed.solve(y=dataNHpi$percDryWgt, Z=msZ, K=aMat, X=msXdb, SE=T)

h2covar <- c(heritability(msOutWWP), heritability(msOutDWPM), heritability(msOutPDW))
names(h2covar) <- c("wetWgtPlot", "dryWgtPerM", "percDryWgt")

locMeans <- msOutWWP$beta[17:22]
locMeans <- locMeans - min(locMeans)
locMeansSE <- msOutWWP$beta.SE[17:22]

msOutWWPh <- mixed.solve(y=log(dataNHpi$wetWgtPlot+1), Z=msZ, K=hMat, X=msX, SE=T)
msOutDWPMh <- mixed.solve(y=log(dataNHpi$dryWgtPerM+1), Z=msZ, K=hMat, X=msX, SE=T)
msOutPDWh <- mixed.solve(y=dataNHpi$percDryWgt, Z=msZ, K=hMat, X=msX, SE=T)

h2covarH <- c(heritability(msOutWWPh), heritability(msOutDWPMh), heritability(msOutPDWh))
names(h2covarH) <- c("wetWgtPlot", "dryWgtPerM", "percDryWgt")

locMeansH <- msOutWWPh$beta[17:22]
locMeansH <- locMeansH - min(locMeansH)
locMeansSEh <- msOutWWPh$beta.SE[17:22]

##########################
# Inidividual measurements
# Somewhere here I would have to filter dataNHim to retain only the ten with
# longest blades
dataNHim$plotSub <- paste(dataNHim$plotNo, dataNHim$subPlot, sep="_")

for (col in c("bladeLength", "bladeMaxWidth", "bladeWidth10cm", "bladeThickness", "stipeLength", "stipeDiameter")) dataNHim[,col] <- as.numeric(dataNHim[,col])

for (col in c("plotNo", "subPlot", "plotSub")) dataNHim[,col] <- as.factor(dataNHim[,col])

# Complicated by the fact that there are no data for some im plots
emZpl <- model.matrix( ~ -1 + plotNo, data=dataNHim)
levels(dataNHim$plotNo) <- levels(dataNHpi$plotNo)
msXim <- model.matrix( ~ -1 + plotNo, data=dataNHim)
msZim <- emZsp <- msXim %*% msZ
msXim <- msXim %*% msXdb

# bladeLength, bladeMaxWidth, bladeWidth10cm, bladeThickness, stipeLength, stipeDiameter, stipeHollow, bullations, fouling
msOutBLh <- mixed.solve(y=log(dataNHim$bladeLength+1), Z=msZim, K=hMat, X=msXim, SE=T)
msOutBMWh <- mixed.solve(y=log(dataNHim$bladeMaxWidth+1), Z=msZim, K=hMat, X=msXim, SE=T)
msOutBW10h <- mixed.solve(y=log(dataNHim$bladeWidth10cm+1), Z=msZim, K=hMat, X=msXim, SE=T)
msOutBTh <- mixed.solve(y=log(dataNHim$bladeThickness+1), Z=msZim, K=hMat, X=msXim, SE=T)
msOutSLh <- mixed.solve(y=log(dataNHim$stipeLength+1), Z=msZim, K=hMat, X=msXim, SE=T)
msOutSDh <- mixed.solve(y=log(dataNHim$stipeDiameter+1), Z=msZim, K=hMat, X=msXim, SE=T)
h2covarHim <- c(heritability(msOutBLh), heritability(msOutBMWh), heritability(msOutBW10h), heritability(msOutBTh), heritability(msOutSLh), heritability(msOutSDh))
names(h2covarHim) <- c("bladeLength", "bladeMaxWidth", "bladeWidth10cm", "bladeThickness", "stipeLength", "stipeDiameter")

# Can't remember exactly why I wanted to do a multikernel.  I think having
# to do with the different levels on individuals and plots.  But I think
# the model did not converge
emZlist <- list(emZsp, emZpl)
emKlist <- list(hMat, diag(ncol(emZpl)))
emOutBLh <- EMMREML::emmremlMultiKernel(y=log(dataNHim$bladeLength+1), X=msXim, Zlist=emZlist, Klist=emKlist)

allBLUPs <- cbind(msOutWWPh$u, msOutDWPMh$u, msOutPDWh$u, msOutBLh$u, msOutBMWh$u, msOutBW10h$u, msOutBTh$u, msOutSLh$u, msOutSDh$u)
colnames(allBLUPs) <- c("WWP", "DWpM", "PDW", "BLen", "BMax", "B10", "BThk", "SLen", "SDia")
corrplot::corrplot.mixed(cor(allBLUPs), diag="n", tl.cex=0.9, tl.col=1)

fndNames <- rownames(allBLUPs)[1:195]
hasDesc <- sapply(fndNames, function(n) grep(n, rownames(allBLUPs)[311:509]))
nDesc <- sapply(hasDesc, length)
colSPGP <- c(rep("black", 195), rep("dark green", 115), rep("dark red", 199))
colSPGP[which(nDesc == 0)] <- "grey"
plot(allBLUPs[,2], pch=16, col=colSPGP, xlab="Individual number", ylab="Dry weight per plot")
rnPos <- rownames(allBLUPs)[intersect(which(nDesc == 0), which(allBLUPs[,2]>0))]
rnNeg <- rownames(allBLUPs)[intersect(which(nDesc == 0), which(allBLUPs[,2]<0))]
locPos <- table(substring(rnPos, 6, 7))
locNeg <- table(substring(rnNeg, 6, 7))

allBLUPsDF <- as.data.frame(allBLUPs)
allBLUPsDF <- cbind(pedigree=rownames(allBLUPs), allBLUPsDF)
rownames(allBLUPsDF)[311:509] <- paste("SL18-UCONN-S", dataNHpi$plotNo[1:199], sep="") 

pdf("DWpMvsPDWsporophytes.pdf")
plot(allBLUPsDF$PDW[311:509], allBLUPsDF$DWpM[311:509], pch=16, cex=0.8, ylab="Dry weight per meter", xlab="Percent dry weight", main="Sporophytes", cex.lab=1.3)
dev.off()
pdf("DWpMvsPDWgametophytes.pdf")
plot(allBLUPsDF$PDW[196:310], allBLUPsDF$DWpM[196:310], pch=16, cex=0.8, ylab="Dry weight per meter", xlab="Percent dry weight", main="Gametophytes", cex.lab=1.3)
dev.off()

# Add gametophyte sex as a column to the dataframe
getGPsex <- function(gpName){
  sexNum <- strsplit(gpName, "-")[[1]][4]
  sex <- substring(sexNum, 1, 1)
  return(ifelse(sex %in% c("F", "M"), sex, NA))
}
gpSex <- sapply(rownames(allBLUPsDF), getGPsex)

# A selection index for sporophytes and gametophytes that weights DWpM twice as high as PDW
allBLUPsDF <- cbind(pedigree=allBLUPsDF$pedigree, gpSex=gpSex, index=2*allBLUPsDF$DWpM + allBLUPsDF$PDW, allBLUPsDF[,-1])
saveRDS(allBLUPsDF, file="allBLUPsDF_analyzeNH.rds")

# Best sporophytes
bestSP_DWpM <- allBLUPsDF[311:509,][order(allBLUPsDF$DWpM[311:509], decreasing=T)[1:20],]
bestSP_DWpM[,-(1:2)] <- bestSP_DWpM[,-(1:2)] %>% round(3)
write.table(bestSP_DWpM, "BestSPbyDryWgtPerM.txt", quote=F, row.names=T, col.names=T)
bestSP_idx <- allBLUPsDF[311:509,][order(allBLUPsDF$index[311:509], decreasing=T)[1:20],]
bestSP_idx[,-(1:2)] <- bestSP_idx[,-(1:2)] %>% round(3)
write.table(bestSP_idx, "BestSPbyDWpMandPDW.txt", quote=F, row.names=T, col.names=T)

# Best gametophytes
bestGP_DWpM <- allBLUPsDF[196:310,][order(allBLUPsDF$DWpM[196:310], decreasing=T)[1:20],]
bestGP_DWpM[,-(1:2)] <- bestGP_DWpM[,-(1:2)] %>% round(3)
write.table(bestGP_DWpM[,-1], "BestGPbyDryWgtPerM.txt", quote=F, row.names=T, col.names=T)
bestGP_idx <- allBLUPsDF[196:310,][order(allBLUPsDF$index[196:310], decreasing=T)[1:20],]
bestGP_idx[,-(1:2)] <- bestGP_idx[,-(1:2)] %>% round(3)
write.table(bestGP_idx[,-1], "BestGPbyDWpMandPDW.txt", quote=F, row.names=T, col.names=T)

# Ordered female and male gametophytes
temp <- allBLUPsDF[allBLUPsDF$gpSex == "F",-(1:2)]
temp <- temp[order(temp$DWpM, decreasing=T),] %>% round(3)
write.table(temp, "FemaleGP_OrderedByDryWgtPerM.txt", quote=F, row.names=T, col.names=T)

temp <- allBLUPsDF[allBLUPsDF$gpSex == "M",-(1:2)]
temp <- temp[order(temp$DWpM, decreasing=T),] %>% round(3)
write.table(temp, "MaleGP_OrderedByDryWgtPerM.txt", quote=F, row.names=T, col.names=T)
