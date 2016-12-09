#' checkCuts
#'
#' Determines the sequence around a cut site using a fasta file or BSgenome
#'
#' @param cutSites A GRanges object containing the locations of the cut sites to be checked for sequence match.
#' @param cutIDs A character vector containing unique IDs for each cut site.
#' @param genome The path to a fasta file or a BSgenome object to check for genomic sequences.
#' @param fasta TRUE if a fasta file has been supplied. Default = FALSE
#' @param seq The desired recognition sequence that the enzyme should have cut.
#' @usage checkCuts(cutSites, cutIDs, genome, fasta = FALSE, seq)
#' @author Benjamin Mayne
#' @return cutIDs that match the desired recognition sequence (seq)
#' @importFrom R.utils gunzip gzip
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @examples
#' library(GenomicRanges)
#' # Load the positions of possible MspI cut sites
#' data(datCounts)
#' cutSites <- GRanges(rownames(datCounts))
#' # Adjust positions to overlap the recognition sequence of MspI
#' start(cutSites) <- start(cutSites) - 1
#' end(cutSites) <- end(cutSites) + 2
#' chr20 <- system.file("extdata", "chr20.fa.gz", package="msgbsR")
#' correctCuts <- checkCuts(cutSites = cutSites, cutIDs = row.names(datCounts),
#'                          genome = chr20, fasta = TRUE, seq = "CCGG")
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

  if(isTRUE(fasta)){
    ## Check if the fasta file is compressed or not
    extenstion <- function (x)
    {
      pos <- regexpr("\\.([[:alnum:]]+)$", x)
      ifelse(pos > -1L, substring(x, pos + 1L), "")
    }
    fastaExt <- ifelse(extenstion(genome) == 'gz', yes=TRUE, no=FALSE)
    if(isTRUE(fastaExt)){
      print('Uncompressing fasta file')
      gunzip(genome)
      genome <- gsub('.gz', '', genome)
    }
  }

  if(!isTRUE(fasta)){
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
  if(isTRUE(fasta)){
    # Scan the fasta file for the sequences
    sequences <- data.frame(scanFa(genome, cutSites))

    if(fastaExt == 'TRUE'){
      print('Compressing fasta file')
      gzip(genome)
    }


  } else {
    sequences <- data.frame(getSeq(genome, cutSites))
  }

  # Add the cutIDs into the data frame
  sequences$ID <- cutIDs

  # Subset for cut sites with the desired sequence
  sequences <- sequences[which(sequences[,1] == seq),]

  # Return the correct sequences
  return(sequences[,2])

}
