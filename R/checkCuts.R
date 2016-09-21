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

  #=====================================================================================

  ## The get.fasta function from the bedr package is below and it has been modified to turn verbose = FALSE
  ## such that there are no messages

  getFASTA <- function (x, fasta = NULL, bed12 = FALSE, strand = FALSE, output.fasta = FALSE,
                        use.name.field = FALSE, check.zero.based = TRUE, check.chr = TRUE,
                        check.valid = TRUE, check.sort = TRUE, check.merge = TRUE,
                        verbose = FALSE)
  {
    catv("FASTA-QUERY\n")
    if (is.null(fasta)) {
      fasta <- system.file("extdata/ucsc.hg19.example.fasta",
                           package = "bedr")
    }
    if (!file.exists(fasta)) {
      catv(" * Fasta does not exist... FAIL\n")
    }
    if (!file.exists(paste0(fasta, ".fai"))) {
      catv(" * Fasta index '.fai' does not exist... FAIL\n")
    }
    params <- paste("-fi", fasta)
    if (use.name.field) {
      params <- paste(params, "-name")
    }
    if (strand) {
      params <- paste(params, "-s")
    }
    if (bed12) {
      params <- paste(params, "-split")
    }
    if (!output.fasta) {
      params <- paste(params, "-tab")
    }
    params <- paste(params, "-fo stdout")
    if (check.valid) {
      is.valid <- is.valid.region(x, check.zero.based = check.zero.based,
                                  check.chr = check.chr, throw.error = TRUE, verbose = verbose)
    }
    region.size <- size.region(x, check.zero.based = FALSE, check.chr = FALSE,
                               check.valid = FALSE, verbose = FALSE)
    if (any(region.size > 8095)) {
      outputFile <- tempfile(pattern = paste("get.fasta", sep = ""),
                             fileext = ".bed")
    }
    else {
      outputFile <- NULL
    }
    x <- bedr(engine = "bedtools", input = list(bed = x), method = "getfasta",
              outputFile = outputFile, params = params, check.zero.based = FALSE,
              check.chr = FALSE, check.valid = FALSE, check.sort = check.sort,
              check.merge = check.merge, verbose = FALSE)
    if (!output.fasta) {
      colnames(x)[1:2] <- c("index", "sequence")
    }
    return(x)
  }


  size.region <- function(x, zero.based = TRUE, check.zero.based = TRUE, check.chr = TRUE, check.valid = TRUE, verbose = TRUE) {

    x <- convert2bed(x, check.zero.based = check.zero.based, check.chr = check.chr, check.valid = check.valid, check.sort = FALSE, check.merge = FALSE, verbose = verbose);

    size <- x$end - x$start;
    if (!zero.based) size <- size + 1;

    return(size)
  }

  #============================================================================================

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

  ## Check if the fasta file is compressed or not
    extenstion <- function (x)
  {
    pos <- regexpr("\\.([[:alnum:]]+)$", x)
    ifelse(pos > -1L, substring(x, pos + 1L), "")
  }
  fastaExt <- ifelse(extenstion(chr1) == 'gz', yes=TRUE, no=FALSE)
  if(fastaExt == 'TRUE'){
    stop('fasta file is compressed, please uncompress the file before using checkCuts')
  }

  ## Check if the seq is a character class
  seqClass <- class(seq)[1]
  seqClass <- ifelse(seqClass == 'character', yes=TRUE, no=FALSE)
  if(seqClass == 'FALSE'){
    stop('seq must be of character class')
  }

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

  print('sorting the cut sites')
  # Sort the cutSites
  if(chr.prefix == TRUE){
    cutSites_sorted <- bedr.sort.region(cutSites, verbose = FALSE)
  } else {
    cutSites_sorted <- bedr.sort.region(cutSites, check.chr = FALSE, verbose = FALSE)
  }

  print('running fasta query')
  # Use the cutSites_sorted to get the sequence of each cut site along with the supplied fasta file
  if(chr.prefix == TRUE){
    sequences <- getFASTA(cutSites_sorted, fasta=fastaPath,
                                         verbose = FALSE)
  } else {
    sequences <- getFASTA(cutSites_sorted, fasta=fastaPath,
                                                check.chr = FALSE, verbose = FALSE)
  }


  # Filter out cutSites that do not match the seq given
  sequences <- sequences[which(sequences[ ,2] == seq | sequences[ ,2] == toupper(seq)), ]

  # Stop function if none of the positions match the sequence
  seqRows <- ifelse(nrow(sequences) == 0, yes = TRUE, no=FALSE)
  if(seqRows == 'TRUE'){
    return('0 cut sites matched the sequence')
  } else {

  row.names(sequences) <- NULL

  # Match the cutSites with thier IDs
  sequences$ID <- cutIDs$cutIDs[match(sequences[,1], cutIDs[,2])]
  sequences <- sequences[,c(3,1,2)]

  # Return the correct sequences
  return(sequences)

  }

}
