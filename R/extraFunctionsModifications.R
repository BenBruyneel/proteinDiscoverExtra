library(stringr)
library(purrr)

#' Helper function for the location data in the different tables. Gets the start
#'  location from a string with a format like: "[100-400]"
#' 
#' @param theString a character vector with the above mentioned format
#'
#' @return integer number (the number before -)
#' @noRd
giveFirst <- function(theString){
  posmin <- stringr::str_locate(theString,pattern = "-")[1,1]
  return(as.integer(substr(theString,2,posmin-1)))
}

#' Helper function for the location data in the different tables. Gets the end
#'  location from a string with a format like: "[100-400]"
#' 
#' @param theString a character vector with the above mentioned format
#'
#' @return integer number (the number after -)
#' @noRd
giveLast <- function(theString){
  posmin <- stringr::str_locate(theString,pattern = "-")[1,1]
  return(as.integer(substr(theString,posmin+1,nchar(theString)-1)))
}

#' function that 'translates' a modification string into a data.frame with three
#'  columns: aa (amino acid), position and modification. If the position is NA,
#'  then that info is not available. This function ONLY works on the
#'  Modification field in the PSMS table. This is because there are differences
#'  between how modificationss are 'encoded' in the other tables
#'
#' @param theString single element vector with info on present modifications
#' @param aaLeaveOut if TRUE then the amino acid info is left out of the
#'  resulting data.frame
#' @param splitchars defines how the modifications are separated in the 
#'  modification string
#' @param patternchars defines the beginning character of a modification name
#'
#' @note this function expects the format of the modification string to be
#'  something like "K5(modification1); C65(modification2); N-Term(Prot)(Acetyl)" 
#' @note it's very well possible that there are modification strings that this
#'  function will have trouble with. It works good for everything the author has
#'  encountered so far. Obviously it will be adjusted whne needed
#'
#' @return a data.frame
#' @export
mods <-function(theString, aaLeaveOut = FALSE, splitchars = "; ", patternchars = "\\("){
  theString <- unlist(strsplit(theString, split=splitchars))
  theString <- lapply(theString,
                      function(x){
                        if (x == "" | purrr::is_empty(x)){
                          return(data.frame(aa = as.character(),
                                            position = as.integer(),
                                            modification = as.character()))
                        }
                        tempi <- stringr::str_locate_all(x,
                                                         pattern = patternchars)
                        if (!is.list(tempi)){
                          tempi <- list(tempi)
                        }
                        tempi <- unname(tempi[[1]][nrow(tempi[[1]]),1])
                        modi <- substr(x, tempi+1, nchar(x)-1)
                        if ((substr(x,1,tempi-1)) == "N-Term(Prot)"){
                          posi <- "1"
                          aai <-"-"
                        } else {
                          posi <- substr(x, 2,tempi-1)
                          aai <- substr(x,1,1)
                        }
                        if (aaLeaveOut){
                          return(data.frame(position = posi,
                                            modification = modi))
                        } else {
                          return(data.frame(aa = aai,
                                            position = posi,
                                            modification = modi))
                        }
                      }
  )
  if (length(theString) == 0){
    return(data.frame(aa = as.character(),
                      position = as.integer(),
                      modification = as.character()))
  }
  tempdf <- theString[[1]]
  if (length(theString) > 1){
    for (counter in 2:((length(theString)))){
      tempdf <- bind_rows(tempdf, theString[[counter]])
    }
  }
  if (class(tempdf$position) != "numeric"){
    for (counter in 1:nrow(tempdf)){
      if ((tempdf$aa[counter] == "N") & (grepl(tempdf$position[counter], pattern = "Term"))){
        tempdf$position[counter] <- 1
      } else {
        # for dealing with peptide C-term mods, untested
        # will need to set all peptides with peptideLocation == -1 to
        # peptideLocation = nchar(peptide)
        # for future expansion of function(s)
        if ((tempdf$aa[counter] == "C") & (grepl(tempdf$position[counter], pattern = "Term"))){
          tempdf$position[counter] <- -1
        }
      }
    }
  }
  return(tempdf %>% mutate(position = as.integer(position)) %>% arrange(position, modification))
}


#' function that combines a set of psms style modifications in either a
#'  set of modification combination strings present or a single string
#'  describing all possible combinations. Essentially this allows to generate
#'  a character vector that describes all modifications present in a peptide
#'  (source of all PSMS's). The combinations/modifcations are only the ones
#'  which have 'proof' in the PSMS table
#'
#' @param theStrings the psms modifications as a character vector
#' @param splitchars defines how the modifications are separated in the 
#'  modification string
#' @param patternchars defines the beginnin
#' @param singleString if TRUE than a single element character vector of all
#'  found modifications is the output, otherwise the different modification
#'  (combinations) will be separate elements
#'
#' @return character vector
#' @export
combinedModsToString <- function(theStrings, splitchars = "; ", patternchars = "\\(", singleString = FALSE){
  if (length(theStrings) < 2){
    return(theStrings)
  }
  modsList <- lapply(theStrings,function(x){mods(x,splitchars = splitchars, patternchars = patternchars)})
  tempdf <- full_join(modsList[[1]], modsList[[2]], by = "position")
  tempdf <- tempdf %>%
    mutate(aa = map2_chr(tempdf$aa.x, tempdf$aa.y, ~ifelse(is.na(.x),.y,.x))) %>%
    select(-c("aa.x","aa.y")) %>%
    select(aa, everything())
  if (length(modsList)>2){
    for (counter in 3:(length(modsList))){
      tempdf <- full_join(tempdf, modsList[[counter]], by = "position")
      tempdf <- tempdf %>%
        mutate(aa = map2_chr(tempdf$aa.x, tempdf$aa.y, ~ifelse(is.na(.x),.y,.x))) %>%
        select(-c("aa.x","aa.y")) %>%
        select(aa, everything())
    }
  }
  tempdf[is.na(tempdf)] <- "Unmodified"
  tempdf <- tempdf %>% arrange(position) %>% t() %>% as.data.frame()
  tempStr <- as.character()
  for (counter in 1:(ncol(tempdf))){
    tempStr <- append(tempStr,paste(c(str_trim(tempdf[1,counter], side = "both"),
                                      str_trim(tempdf[2,counter], side = "both"),
                                      "(",paste(sort(unique(tempdf[3:nrow(tempdf),counter])), collapse = "/"),")"), collapse = ""))
  }
  if (singleString){
    return(paste(tempStr, collapse = "; "))
  } else {
    return(tempStr)
  }
}

#' function to figure out if a peptide (identified by PeptideGroupIDs) contains
#'  a modification. Can be done with and without considering the position of 
#'  the modification
#'
#' @param db database access 'handle'
#' @param PeptideIDs the peptideID for which the data is to be retrieved.
#'  This should be a numeric or character vector format. The output from the
#'  dbGetPeptideIDs function gives a data.frame so is only compatible when
#'  using the column "TargetPeptideGroupsPeptideGroupID". Note: this function is
#'  not vectorized.
#' @param modificationName name of the modification. Note that this should be
#'  the modification name as used by the psms table (TargetOsms), not the
#'  modification name used by the modification table (ModificationSites)
#' @param modificationPosition position of the modification in the peptide, if
#'  the modification is present, but not in the modificationPosition, then the
#'  function will return FALSE. This setting is ignored if it is NA
#'
#' @return logical vector (TRUE or FALSE)
#' @export
peptideContainsModification <- function(db, PeptideGroupIDs,
                                        modificationName,
                                        modificationPosition = NA){
  if (identical(modificationName,NA)){
    stop("Error : Modification Name cannot be NA")
  }
  tempMods <- (dbGetPsmIDs(db = db, PeptideGroupIDs = PeptideGroupIDs) %>%
                 dbGetPsmTable(db = db) %>%
                 filter(!is.na(MasterProteinAccessions))
  )$Modifications
  tempMods <- mods(tempMods) %>% distinct(.keep_all = TRUE)
  if (!identical(modificationPosition,NA)){
    tempMods <- tempMods %>% filter(position %in% modificationPosition, modification %in% modificationName)
    return(nrow(tempMods) != 0)
  } else {
    tempMods <- tempMods %>% filter(modification %in% modificationName)
    return(nrow(tempMods) != 0)
  }
}


#' funnction to extract the start and end positions of a peptide in a
#'  specific protein
#'
#' @param peptideTableRow a row from a peptide table containing the
#'  'PositionsinMasterProteins' field. This field may contain more than one
#'   master protein and therefore the function needs an accession number to
#'   specify from which one the positions need to be extracted
#' @param Accession accession id (Uniprot style) of the (master) protein for
#'  which the positions need to be calculated
#' @param multipleAccession integer vector, default = NA. In case of the same
#'  accession occurring more than once, specifies which one to use. This can
#'  happen when multiple databases with overlapping accessions are used
#' @param showWarning logical vector: give warning when multiple accession
#'  occurs
#' @param removeNA specifies what to to do when the result is NA: if TRUE
#'  returns an empty peptideTableRow, if FALSE then returns the original
#'  peptideTableRow with NA in the two extra columns (starLocation &
#'  endLocation)
#' @param positionsIn specifies the column to get the positions from. This can
#'  be different depending on Proteome Discoverer settings. Can also be
#'  'PositionsinProteins'
#'
#' @return a data.frame row with the same data as the argument peptideTableRow
#'  plus the start & end position of the peptide (present)
#' @export
givePositions <- function(peptideTableRow, Accession,
                          multipleAccession = NA, showWarning = FALSE,
                          removeNA = FALSE,
                          positionsIn = "PositionsinMasterProteins"){
  # remove all [  ] characters
  #tempStr <- str_replace_all(str_replace_all(peptideTableRow$PositionsinMasterProteins, pattern = "\\[", replacement = ""), pattern = "\\]", replacement = "")
  tempStr <- peptideTableRow[, positionsIn]
  tempStr <- strsplit(tempStr, split = "; ")[[1]]
  tempStr <- str_trim(tempStr, side = "both")
  # find strings that have accessions 
  whichStr <- which(str_detect(tempStr, pattern = " \\[\\d+-\\d+\\]"))
  # find which string has the accession
  thisIsAccession <- whichStr[which(grepl(tempStr[whichStr], pattern = Accession))]
  # end of Accession
  if (!is_empty(thisIsAccession)){
    # in case the accession exists more than one, possible with multiple database searches!
    if (length(thisIsAccession)>1){
      if (!is.na(multipleAccession)){
        thisIsAccession <- 1       # thisIsAccession[multipleAccession] # remove
        tempStr <- tempStr[multipleAccession]
      }
      if (showWarning){
        warning(" Multiple identical accessions in peptide table row!")
      }
    }
    endOfAccession <- whichStr[which(grepl(tempStr[whichStr], pattern = Accession))+1]-1
    if (!identical(endOfAccession, NA)){
      endOfAccession <- length(tempStr)
    }
    tempStr <- tempStr[thisIsAccession:endOfAccession]
    tempdf <- data.frame()
    for (counter in 1:length(tempStr)){
      newdf <- peptideTableRow
      tempStr[counter] <- unlist(str_extract(tempStr[counter],
                                             pattern = "\\[\\d+-\\d+\\]"))
      newdf$startLocation <- giveFirst(tempStr[counter])
      newdf$endLocation <- giveLast(tempStr[counter])
      tempdf <- bind_rows(tempdf, newdf)
    }
    return(tempdf)
  } else {
    peptideTableRow$startLocation <- NA
    peptideTableRow$endLocation <- NA
    if (removeNA){
      return(peptideTableRow[0,])
    } else {
      return(peptideTableRow)
    }
  }
}


#' function to either calculate the percentage of a certain modification present
#'  at a specific protein position or to generate a table of data on all
#'  peptides in the protein containing the specific protein position
#'
#' @param db database access 'handle'
#' @param peptideTable peptide table coming from eg
#'  \code{\link{dbGetPeptideTable}} 
#' @param modificationLocation position in protein (not peptide!) for which
#'  to get or calculate data
#' @param modificationName name of the modification, can be not present at
#'  position. Note that, depending on the modification definition in
#'  proteome discoverer, this argument MUST be the name as used in the psms
#'  table, which is the Abbreviation and NOT the Name of the modification.
#'  The information regarding Name/Abbreviation can be found in the
#'  'FoundModifications' table in a pdResult file/database
#' @param Accession accession id (Uniprot style) of the (master) protein
#' @param multipleAccession integer vector, default = NA. In case of the same
#'  accession occurring more than once, specifies which one to use. This can
#'  happen when multiple databases with overlapping accessions are used
#' @param showWarning logical vector: give warning when multiple accession
#'  occurs
#' @param positionsIn specifies the column to get the positions from. This can
#'  be different depending on Proteome Discoverer settings. Can also be
#'  'PositionsinProteins'
#' @param giveTable if FALSE then the result will be a single row data.frame
#'  with the modificationName (Abbreviation) of the modification, the location
#'  in protein and the percentage of label (based on abundance of all peptides
#'  containing the location). If TRUE then the result will be a data.frame with
#'  all the 'raw' information. It will contain all peptides which contain the
#'  position (with and without modification), abundances, etc etc. It will also
#'  have a column specifying whether a row/peptide contains the modification
#'
#' @return data.frame
#' @export
calculatePositionPercentage <- function(db, peptideTable, modificationLocation,
                                        modificationName, Accession,
                                        multipleAccession = NA, showWarning = FALSE,
                                        positionsIn = "PositionsinMasterProteins",
                                        giveTable = FALSE){
  tempdf <- data.frame()
  for (counter in 1:nrow(peptideTable)){
    newdf <- givePositions(peptideTableRow = peptideTable[counter,],
                           Accession = Accession,
                           multipleAccession = multipleAccession,
                           showWarning = showWarning,
                           removeNA = TRUE,
                           positionsIn = positionsIn)
    tempdf <- bind_rows(tempdf, newdf)
  }
  peptideTable <- tempdf
  exa <- peptideTable %>%
    filter(startLocation <= modificationLocation &
             endLocation >= modificationLocation) %>%
    select(Sequence, Modificationsallpossiblesites, Abundances_1,
           PeptideGroupID, startLocation, endLocation)
  if (nrow(exa)>0){
    exa$peptideLocation <- modificationLocation - exa$startLocation + 1
    exa$modificationPresent <- FALSE
    for (counter in 1:nrow(exa)){
      exa$modificationPresent[counter] <- 
        peptideContainsModification(db = db,
                                    PeptideGroupIDs = exa$PeptideGroupID[counter],
                                    modificationName = modificationName,
                                    modificationPosition = exa$peptideLocation[counter])
    }
    if (giveTable){
      return(exa)
    } else {
      return(data.frame(modification = modificationName,
                        location = modificationLocation,
                        percentage = (sum(exa$Abundances_1[exa$modificationPresent],na.rm = TRUE)/
                                        sum(exa$Abundances_1,na.rm = TRUE)) * 100 ))
    }
  } else {
    if (giveTable){
      return(NA)
    } else {
      data.frame(modification = as.character(), location = as.integer(),
                 percentage = as.numeric())
    }
  }
}


#' function to generate a table of data on all peptides in the protein
#'  containing q specific protein position
#'
#' @param db database access 'handle'
#' @param peptideTable peptide table coming from eg
#'  \code{\link{dbGetPeptideTable}} 
#' @param modificationLocation position in protein (not peptide!) for which
#'  to get or calculate data
#' @param modificationName name of the modification, can be not present at
#'  position. Note that, depending on the modification definition in
#'  proteome discoverer, this argument MUST be the name as used in the psms
#'  table, which is the Abbreviation and NOT the Name of the modification.
#'  The information regarding Name/Abbreviation can be found in the
#'  'FoundModifications' table in a pdResult file/database
#' @param Accession accession id (Uniprot style) of the (master) protein
#' @param multipleAccession integer vector, default = NA. In case of the same
#'  accession occurring more than once, specifies which one to use. This can
#'  happen when multiple databases with overlapping accessions are used
#' @param showWarning logical vector: give warning when multiple accession
#'  occurs

#' @param positionsIn specifies the column to get the positions from. This can
#'  be different depending on Proteome Discoverer settings. Can also be
#'  'PositionsinProteins'
#'
#' @note Essentially this function does the same as
#'  \code{\link{calculatePositionPercentage}} but does not only calculate
#'  the percentage of a position
#' @note This function is kept for compatibility with other older code and will
#'  likely be removed eventually
#'
#' @return data.frame
#' @export
getPositionTable <- function(db, peptideTable, modificationLocation,
                             modificationName, Accession,
                             multipleAccession = NA, showWarning = FALSE,
                             positionsIn = "PositionsinMasterProteins"){
  tempdf <- data.frame()
  for (counter in 1:nrow(peptideTable)){
    newdf <- givePositions(peptideTableRow = peptideTable[counter,],
                           Accession = Accession,
                           multipleAccession = multipleAccession,
                           showWarning = showWarning,
                           removeNA = TRUE,
                           positionsIn = positionsIn)
    tempdf <- bind_rows(tempdf, newdf)
  }
  peptideTable <- tempdf
  exa <- peptideTable %>%
    filter(startLocation <= modificationLocation &
             endLocation >= modificationLocation) %>%
    select(Sequence, Modificationsallpossiblesites, Abundances_1,
           PeptideGroupID, startLocation, endLocation)
  exa$peptideLocation <- modificationLocation - exa$startLocation + 1
  exa$modificationPresent <- FALSE
  for (counter in 1:nrow(exa)){
    exa$modificationPresent[counter] <- 
      peptideContainsModification(db = db,
                                  PeptideGroupIDs = exa$PeptideGroupID[counter],
                                  modificationName = modificationName,
                                  modificationPosition = exa$peptideLocation[counter])
  }
  return(exa)
}


#' Function that takes the peptide table info from
#'  \code{\link{getPositionTable}} or \code{\link{calculatePositionPercentage}}
#'  and calculates the precentage of modified peptides (specific modification,
#'  specific location)
#'
#' @param positionTable result form either \code{\link{getPositionTable}} or
#'  \code{\link{calculatePositionPercentage}}. Must have 'Abundance', location
#'  and modificationPresent information
#' @param modificationLocation position in protein (not peptide!) for which
#'  to get or calculate data
#' @param modificationName name of the modification. This argument is in this
#'  function only used in the resulting data.frame
#'
#' @return data.frame (modification, location, percentage)
#' @export
getPositionPercentage <- function(positionTable,
                                  modificationLocation, modificationName){
  return(data.frame(modification = modificationName,
                    location = modificationLocation,
                    percentage = (sum(positionTable$Abundances_1[positionTable$modificationPresent],
                                      na.rm = TRUE) /
                                    sum(positionTable$Abundances_1,na.rm = TRUE)) * 100 ))
}

#' function that generates a complete 
#'
#' @param db database access 'handle'
#' @param modificationTable data.frame coming from
#'  \code{\link{dbGetModificationTable}}. Obviously the protein for which the
#'  modification was retrieved MUST match the Accession argument. Also: it is
#'  necessary to replace the modification Name with the Abbreviation as used
#'  in the FoundModifications table (if the Name and the Abbreviation are not
#'  the same!)
#' @param peptideTable peptide table coming from eg
#'  \code{\link{dbGetPeptideTable}} 
#' @param Accession accession id (Uniprot style) of the (master) protein for
#'  which the table is to be generated
#' @param multipleAccession integer vector, default = NA. In case of the same
#'  accession occurring more than once, specifies which one to use. This can
#'  happen when multiple databases with overlapping accessions are used
#' @param showWarning logical vector: give warning when multiple accession
#'  occurs
#' @param positionsIn specifies the column to get the positions from. This can
#'  be different depending on Proteome Discoverer settings. Can also be
#'  'PositionsinProteins'
#'
#' @return data.frame of all modifications, their locations in the protein and
#'  the percentage of each per position
#' @export
proteinModificationTable <- function(db, modificationTable, peptideTable,
                                     Accession,
                                     multipleAccession = NA, showWarning = FALSE,
                                     positionsIn = "PositionsinMasterProteins"){
  modsperc <- map2_df(modificationTable$ModificationName,
                      modificationTable$Position,
                      ~calculatePositionPercentage(db = db,
                                                   peptideTable = peptideTable,
                                                   modificationLocation = .y,
                                                   modificationName = .x,
                                                   Accession = Accession,
                                                   multipleAccession = multipleAccession,
                                                   showWarning = showWarning,
                                                   positionsIn = positionsIn))
  return(modsperc)
}

#' @Description function around stringr::str_locate_all function to retrieve ONLY the start locations of a vector
#'
#' @note not vectorized, for more info on the parameters see ?stringr::str_locate_all
#' @param string Input vector. Either a character vector, or something coercible to one
#' @param pattern pattern to look for
#'
#' @return numeric vector
#' @export
str_locate_all_starts <- function(string, pattern){
  return(str_locate_all(string = string, pattern = pattern)[[1]][,"start"])
}

#' @Description function around stringr::str_locate_all function to retrieve ONLY the end locations of a vector
#'
#' @note not vectorized, for more info on the parameters see ?stringr::str_locate_all
#' @param string Input vector. Either a character vector, or something coercible to one
#' @param pattern pattern to look for
#'
#' @return numeric vector
#' @export
str_locate_all_ends <- function(string, pattern){
  return(str_locate_all(string = string, pattern = pattern)[[1]][,"end"])
}
