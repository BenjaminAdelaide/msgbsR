#' checkCuts
#'
#' Determines the sequence around a cut site
#'
#' @param cutSites A matrix, where the first column is the chromosome ID, the second column is the starting position and the third column is the ending position.
#' @param fastaPath The path to fasta file to be checked against.
#' @param seq The desired sequence that the enzyme should have cut.
#' @usage checkCuts(cutSites = my_cuts, fastaPath = 'path-to-fasta-file', seq = sequence)
#' @examples
#' \dontrun{
#' # Firstly generate a matrix containing some MspI cut sites (sequence=CCGG)
#' set.seed(1)
#' x <- matrix(data=c(rep('chr1', 3), '29736', '120552', '297122',
#'                                    '29736', '120552', '297122'))
#' x[,3] <- as.numeric(x[,2]) + 2
#' x[,2] <- as.numeric(x[,2]) - 2
#' my_cuts <- checkCuts(cutSites=x, fastaPath='/hg38/genome.fa', seq='CCGG')
#' }
#' @author Benjamin Mayne
#' @export

checkCuts <- function(cutSites, fastaPath, seq){

  # Firstly turn the cutSites into the appropiate format to be sorted
  cutSites <- paste(cutSites[ ,1], ':', cutSites[ ,2], '-', cutSites[ ,3], sep='')

  # Sort the cutSites
  cutSites_sorted <- bedr.sort.region(cutSites)

  # Use the cutSites_sorted to get the sequence of each cut site along with the supplied fasta file
  sequences <- get.fasta(cutSites_sorted, fasta=fastaPath)

  # Filter out cutSites that do not match the seq given
  sequences <- sequences[which(sequences[ ,2] == seq | sequences[ ,2] == toupper(seq)), ]
  row.names(sequences) <- NULLs

  # Return the correct sequences
  return(sequences)

}
