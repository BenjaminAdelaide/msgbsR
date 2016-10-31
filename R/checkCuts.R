#' checkCuts
#'
#' Determines the sequence around a cut site using a fasta file
#'
#' @param cutSites A matrix where the first column is the chromosome ID, the second column is the starting position and the third column is the ending position of the cut site sequence.
#' @param cutIDs A character vector containing unique IDs for each cut site.
#' @param fastaPath The path to fasta file to be checked against, if not using a BSgenome.
#' @param BSgenome The name of the species BSgenome to be used, if not using a fasta file.
#' @param seq The desired sequence that the enzyme should have cut.
#' @param chr.prefix Is there a chr prefix in the first column of cutSites. Default is TRUE.
#' @usage checkCuts(cutSites, cutIDs, fastaPath, BSgenome seq)
#' @author Benjamin Mayne
#' @export

checkCuts <- function(cutSites, cutIDs, fastaPath = NULL, BSgenome = NULL, seq){

  # Unit tests
  ## Check if the cutSites are a matrix
  cutSitesClass <- class(cutSites)[1]
  cutSitesClass <- ifelse(cutSitesClass == 'matrix', yes=TRUE, no=FALSE)
  if(cutSitesClass == 'FALSE'){
    stop('cutSites must be a matrix')
  }

  ## Check if the cutIDs are a character class
  cutIDsClass <- class(cutIDs)[1]
  cutIDsClass <- ifelse(cutIDsClass == 'character', yes=TRUE, no=FALSE)
  if(cutIDsClass == 'FALSE'){
    stop('The cutIDs must be of character class')
  }

  if(is.null(fastaPath) == FALSE){
  ## Check if the fasta file is compressed or not
    extenstion <- function (x)
  {
    pos <- regexpr("\\.([[:alnum:]]+)$", x)
    ifelse(pos > -1L, substring(x, pos + 1L), "")
  }
  fastaExt <- ifelse(extenstion(fastaPath) == 'gz', yes=TRUE, no=FALSE)
  if(fastaExt == 'TRUE'){
    print('Uncompressing fasta file')
    gunzip(fastaPath)
    fastaPath <- gsub('.gz', '', fastaPath)
  }
  }

  if(is.null(BSgenome) == FALSE){
  ## Check if supplied BSgenome object is actually a BSgenome
  output = class(BSgenome) == 'BSgenome'
  if(output == FALSE){
    stop('Supplied input is not a BSgenome')
  }
  }

  ## Check if the seq is a character class
  seqClass <- class(seq)[1]
  seqClass <- ifelse(seqClass == 'character', yes=TRUE, no=FALSE)
  if(seqClass == 'FALSE'){
    stop('seq must be of character class')
  }

#=================================================================================

  # Firstly turn the cutSites into a data frame and then into a GRanges object
  cutSites <- data.frame(cutSites)
  colnames(cutSites) <- c('seqname', 'start', 'end')

  # Check if any numbers have been short to e+etc
  es <- grep('e', cutSites[,3])

  if(length(es) > 0){
    cutSites[es,3] <- as.character(lapply(1:length(es), function(x) sprintf("%.f", as.numeric(cutSites[es[x] ,3]))))
  }
  ## Repeat it for the start column
  es <- grep('e', cutSites[,2])

  if(length(es) > 0){
    cutSites[es,2] <- as.character(lapply(1:length(es), function(x) sprintf("%.f", as.numeric(cutSites[es[x] ,2]))))
  }

  # make sure the start and end are numeric
  cutSites$start <- as.numeric(as.character(cutSites$start))
  cutSites$end <- as.numeric(as.character(cutSites$end))

  # Turn into a GRanges object
  cuts.gr <- makeGRangesFromDataFrame(cutSites)

  # Use either a fasta file or BSgenome
  if(is.null(fastaPath) == FALSE){
  # Scan the fasta file for the sequences
  sequences <- data.frame(scanFa(fastaPath, cuts.gr))

  if(fastaExt == 'TRUE'){
    print('Compressing fasta file')
    gzip(fastaPath)
  }


  } else if(is.null(BSgenome) == FALSE) {
  sequences <- data.frame(getSeq(BSgenome, cuts.gr))
  }

  # Add the cutIDs into the data frame
  sequences$ID <- cutIDs

  # Subset for cut sites with the desired sequence
  sequences <- sequences[which(sequences[,1] == seq),]

  # Return the correct sequences
  return(sequences[,2])

  }
