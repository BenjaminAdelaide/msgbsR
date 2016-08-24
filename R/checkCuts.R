#' checkCuts
#'
#' Determines the sequence around a cut site
#'
#' @param cutSites A matrix where the first column is the chromosome ID, the second column is the starting position and the third column is the ending position.
#' @param cutIDs A character vector containing unique IDs for each cut site.
#' @param fastaPath The path to fasta file to be checked against.
#' @param seq The desired sequence that the enzyme should have cut.
#' @param chr.prefix Is there a chr prefix in the first column of cutSites. Default is TRUE.
#' @usage checkCuts(cutSites = my_cuts, cutIDs = my_IDs, fastaPath = 'path-to-fasta-file', seq = sequence, chr.prefix = TRUE)
#' @examples
#' \dontrun{
#' # Firstly generate a matrix containing some MspI cut sites (sequence=CCGG)
#' x <- matrix(data=c(rep('chr1', 3), '29736', '120552', '297122',
#'                                    '29736', '120552', '297122'))
#' x[,3] <- as.numeric(x[,2]) + 2
#' x[,2] <- as.numeric(x[,2]) - 2
#' y <- c('A', 'B', 'C')
#' my_cuts <- checkCuts(cutSites=x, cutIDs=y fastaPath='/hg38/genome.fa', seq='CCGG')
#' }
#' @author Benjamin Mayne
#' @export

checkCuts <- function(cutSites, cutIDs, fastaPath, seq, chr.prefix = TRUE){

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

  # Turn the cutSites into the appropiate format to be sorted
  cutSites <- paste(cutSites[ ,1], ':', cutSites[ ,2], '-', cutSites[ ,3], sep='')

  # Combine the cutSites with the cutIDs to return them at the end of the function
  cutIDs <- data.frame(cbind(cutIDs, cutSites))

  # Sort the cutSites
  if(chr.prefix == TRUE){
    cutSites_sorted <- bedr.sort.region(cutSites)
  } else {
    cutSites_sorted <- bedr.sort.region(cutSites, check.chr = FALSE)
  }

  # Use the cutSites_sorted to get the sequence of each cut site along with the supplied fasta file
  if(chr.prefix == TRUE){
    sequences <- get.fasta(cutSites_sorted, fasta=fastaPath)
  } else {
    sequences <- get.fasta(cutSites_sorted, fasta=fastaPath, check.chr = FALSE)
  }

  # Filter out cutSites that do not match the seq given
  sequences <- sequences[which(sequences[ ,2] == seq | sequences[ ,2] == toupper(seq)), ]
  row.names(sequences) <- NULL

  # Match the cutSites with thier IDs
  sequences$ID <- cutIDs$cutIDs[match(sequences[,1], cutIDs[,2])]
  sequences <- sequences[,c(3,1,2)]

  # Return the correct sequences
  return(sequences)

}
