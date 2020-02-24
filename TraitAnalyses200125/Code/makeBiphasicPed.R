#' Make a three-column pedigree based on a data.frame formated as dataNH
#'
#' @param kelpNameColumns A two-column data.frame with female parents in
#' column 1 and male parents in column 2.  These columns have names in the
#' specific format of dataNH
#'  
#' @return A three-column pedigree matrix with row number and parent ID
#' listed as row numbers 
#'
makeBiphasicPed <- function(kelpNameColumns, founderFirst=NULL){
  founder <- unlist(kelpNameColumns)
  founder <- strsplit(founder, split="-", fixed=T)
  founder <- sapply(founder, function(vec) paste(vec[1:3], collapse="-"))
  founder <- unique(founder)
  founder <- c(founderFirst, setdiff(founder, founderFirst))
  nFounders <- length(founder)
  parFemale <- unique(kelpNameColumns[,1])
  parMale <- unique(kelpNameColumns[,2])
  founderParFema <- sapply(strsplit(parFemale, split="-", fixed=T), function(vec) paste(vec[1:3], collapse="-"))
  founderParMale <- sapply(strsplit(parMale, split="-", fixed=T), function(vec) paste(vec[1:3], collapse="-"))
  fndRow <- 1:nFounders
  names(fndRow) <- founder
  femaRow <- nFounders + 1:length(parFemale)
  names(femaRow) <- parFemale
  maleRow <- nFounders + length(parFemale) + 1:length(parMale)
  names(maleRow) <- parMale
  fndPed <- cbind(fndRow, 0, 0)
  femaPed <- cbind(femaRow, fndRow[founderParFema], NA)
  malePed <- cbind(maleRow, fndRow[founderParMale], NA)
  sporProg <- 1:nrow(kelpNameColumns) + nFounders + length(parFemale) + length(parMale)
  progPed <- cbind(sporProg, femaRow[kelpNameColumns[,1]], maleRow[kelpNameColumns[,2]])
  allPed <- rbind(fndPed, femaPed, malePed, progPed)
  rownames(allPed) <- c(founder, parFemale, parMale, paste(kelpNameColumns[,1], kelpNameColumns[,2], sep="x"))
  return(allPed)
}
