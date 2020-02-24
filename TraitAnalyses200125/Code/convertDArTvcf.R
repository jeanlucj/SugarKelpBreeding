# Convert DArT format VCF
# Input: file name
# Output: -1, 0, 1 format matrix with individuals in rows and marker names in columns
convertDArTvcf <- function(fileName){
  dartVecToDosage <- function(dv){
    slashToDos <- function(slash){
      switch(slash, "0/0" = -1, "0/1" = 0, "1/1" = 1, "./." = NA)
    }
    result <- sapply(dv, slashToDos)
    names(result) <- NULL
    return(result)
  }
  vcf <- file(fileName)
  open(vcf, "r")
  indNames <- scan(vcf, skip=2, nlines=1, what=character())[-(1:9)]
  nInd <- length(indNames)
  dartMat <- matrix(scan(vcf, what=character()), ncol=nInd+9, byrow=T)
  close(vcf)
  mrkNames <- apply(dartMat[,1:2], 1, function(v) paste(v, collapse="."))
  dosages <- apply(dartMat[,-(1:9)], 1, dartVecToDosage)
  rownames(dosages) <- indNames
  colnames(dosages) <- mrkNames
  return(dosages)
}
