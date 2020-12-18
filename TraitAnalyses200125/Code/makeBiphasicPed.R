#' Make a three-column pedigree based on a data.frame formated as dataNH
#'
#' @param kelpNameColumns A two-column data.frame with female parents in
#' column 1 and male parents in column 2.  These columns have names in the
#' specific format of dataNH.
#' Name format: SA18-CB-10-F1
#' Two letters for the species, two digits for the sampling year
#' dash two letters for the sampling location dash
#' Sporophyte sample number dash
#' Single letter male or female, single digit gametophyt number
#'
#' @return A three-column pedigree matrix with row number and parent ID
#' listed as row numbers
#'
makeBiphasicPed <- function(kelpNameColumns, founderFirst=NULL){
    # Sporophyte founders
    founder <- unlist(kelpNameColumns)
    founder <- strsplit(founder, split="-", fixed=T)
    founder <- sapply(founder, function(vec) paste(vec[1:3], collapse="-"))
    founder <- unique(founder)
    founder <- c(founderFirst, setdiff(founder, founderFirst))
    nFounders <- length(founder)
    # Make the pedigree matrix for the founders
    fndRow <- 1:nFounders
    names(fndRow) <- founder
    fndPed <- cbind(fndRow, 0, 0)
    
    # Gametophytes coming from founders
    parFemale <- unique(kelpNameColumns[,1])
    parMale <- unique(kelpNameColumns[,2])
    # The SP parents of the gametophytes
    founderParFema <- sapply(strsplit(parFemale, split="-", fixed=T), function(vec) paste(vec[1:3], collapse="-"))
    founderParMale <- sapply(strsplit(parMale, split="-", fixed=T), function(vec) paste(vec[1:3], collapse="-"))
    # Make the pedigree matrix for the GPs from the founders
    femaRow <- nFounders + 1:length(parFemale)
    names(femaRow) <- parFemale
    maleRow <- nFounders + length(parFemale) + 1:length(parMale)
    names(maleRow) <- parMale
    # GPs only have one parent, so second parent is NA
    femaPed <- cbind(femaRow, fndRow[founderParFema], NA)
    malePed <- cbind(maleRow, fndRow[founderParMale], NA)
    
    # Sporophyte progeny of the gametophytes
    sporProg <- 1:nrow(kelpNameColumns) + nFounders + length(parFemale) + length(parMale)
    progPed <- cbind(sporProg, femaRow[kelpNameColumns[,1]], maleRow[kelpNameColumns[,2]])
    
    # Put all the bits of pedigree matrix together
    allPed <- rbind(fndPed, femaPed, malePed, progPed)
    rownames(allPed) <- c(founder, parFemale, parMale, paste(kelpNameColumns[,1], kelpNameColumns[,2], sep="x"))
    return(allPed)
}
