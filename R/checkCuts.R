#' checkCuts
#'
#' Determines the sequence around a cut site using a fasta file
#'
#' @param cutSites A matrix where the first column is the chromosome ID, the second column is the starting position and the third column is the ending position of the cut site sequence.
#' @param cutIDs A character vector containing unique IDs for each cut site.
#' @param genome The path to fasta file or a BSgenome object to check for genomic sequences.
#' @param fasta TRUE if a fasta file has been supplied. Default = FALSE
#' @param seq The desired sequence that the enzyme should have cut.
#' @usage checkCuts(cutSites, cutIDs, genome, fasta = FALSE, seq)
#' @author Benjamin Mayne
#' @return cut sites that match the input sequence
#' @importFrom R.utils gunzip gzip
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @examples
#' library(GenomicRanges)
#' data(datCounts)
#' cutSites <- GRanges(rownames(datCounts))
#' cutSites@ranges@start <- as.integer(cutSites@ranges@start -1)
#' cutSites@ranges@width <- as.integer(cutSites@ranges@width + 3)
#' chr20 <- system.file("extdata", "chr20.fa.gz", package="msgbsR")
#' correctCuts <- checkCuts(cutSites = cutSites, cutIDs = row.names(datCounts),
#'                          genome = chr20, fasta = TRUE, seq = 'CCGG')
#' @export


checkCuts <- function(cutSites, cutIDs, genome, fasta = FALSE, seq){

  # Unit tests
  ## Check if the cutSites is GRanges object
  if(!is(cutSites, "GRanges")){
    stop('The cutSites must be a GRanges object')
  }

  ## Check if the cutIDs are a character class
  if(!is(cutIDs, "character")){
    stop('The cutIDs must be of character class')
  }

  ## If a fasta file has been supplied check

  if(fasta == TRUE){
  ## Check if the fasta file is compressed or not
    extenstion <- function (x)
  {
    pos <- regexpr("\\.([[:alnum:]]+)$", x)
    ifelse(pos > -1L, substring(x, pos + 1L), "")
  }
  fastaExt <- ifelse(extenstion(genome) == 'gz', yes=TRUE, no=FALSE)
  if(fastaExt == 'TRUE'){
    print('Uncompressing fasta file')
    gunzip(genome)
    genome <- gsub('.gz', '', genome)
  }
  }

  if(fasta == FALSE){
  ## Check if supplied genome object is actually a BSgenome
    if(!is(genome, "BSgenome")){
      stop('genome must be a BSgenome if fasta = FALSE')
    }
  }

  ## Check if the seq is a character class
  if(!is(seq, "character")){
    stop('seq must be of character class')
  }

#=================================================================================

  # Use either a fasta file or BSgenome
  if(fasta == TRUE){
  # Scan the fasta file for the sequences
  sequences <- data.frame(scanFa(genome, cutSites))

  if(fastaExt == 'TRUE'){
    print('Compressing fasta file')
    gzip(genome)
  }


  } else if(fasta == FALSE) {
  sequences <- data.frame(getSeq(genome, cutSites))
  }

  # Add the cutIDs into the data frame
  sequences$ID <- cutIDs

  # Subset for cut sites with the desired sequence
  sequences <- sequences[which(sequences[,1] == seq),]

  # Return the correct sequences
  return(sequences[,2])

}



