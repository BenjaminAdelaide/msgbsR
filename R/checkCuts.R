#' checkCuts
#'
#' Determines the sequence around a cut site using a fasta file or BSgenome
#'
#' @param cutSites A GRanges object containing the locations of the cut sites to be checked for sequence match. The names of the correct cut sites will be returned as a GRanges object.
#' @param genome The path to a fasta file or a BSgenome object to check for genomic sequences.
#' @param fasta TRUE if a fasta file has been supplied. Default = FALSE
#' @param seq The desired recognition sequence that the enzyme should have cut.
#' @usage checkCuts(cutSites, genome, fasta = FALSE, seq)
#' @author Benjamin Mayne
#' @return A GRanges object containing the names of the sites that had the correct sequence.
#' @importFrom R.utils gunzip gzip
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom BSgenome getBSgenome
#' @examples
#' library(GenomicRanges)
#' library(SummarizedExperiment)
#' library(BSgenome.Rnorvegicus.UCSC.rn6)
#' # Load the positions of possible MspI cut sites
#' data(ratdata)
#' # Extract the cut sites
#' cutSites <- rowRanges(ratdata)
#' # Adjust for positive strand
#' start(cutSites[strand(cutSites) == "+"]) <- start(cutSites[strand(cutSites) == "+"]) - 1
#' end(cutSites[strand(cutSites) == "+"]) <- end(cutSites[strand(cutSites) == "+"]) + 2
#' # Adjust for negative strand
#' start(cutSites[strand(cutSites) == "-"]) <- start(cutSites[strand(cutSites) == "-"]) - 2
#' end(cutSites[strand(cutSites) == "-"]) <- end(cutSites[strand(cutSites) == "-"]) + 1
#' correctCuts <- checkCuts(cutSites = cutSites, genome = "rn6", seq = "CCGG")
#' @export

checkCuts <- function(cutSites, genome, fasta = FALSE, seq){

    # Unit tests
    ## Check if the cutSites is GRanges object
    if(!is(cutSites, "GRanges")){
        stop("The cutSites must be a GRanges object")
    }

    ## If a fasta file has been supplied check

    if(fasta){
        ## Check if the fasta file is compressed or not
        extenstion <- function (x)
        {
            pos <- regexpr("\\.([[:alnum:]]+)$", x)
            ifelse(pos > -1L, substring(x, pos + 1L), "")
        }
        if(extenstion(genome) == "gz"){
            fastaExt = TRUE
            print("Uncompressing fasta file")
            gunzip(genome)
            genome <- gsub(".gz", "", genome)
        }
    }

    if(!fasta){
        ## get BSgenome
        genome <- getBSgenome(genome)

    }

    ## Check if the seq is a character class
    if(!is(seq, "character")){
        stop("seq must be of character class")
    }

    #=================================================================================

    # Use either a fasta file or BSgenome
    if(fasta){
        # Scan the fasta file for the sequences
        sequences <- data.frame(scanFa(genome, cutSites))
        # Add the cutIDs into the data frame
        sequences$ID <- names(cutSites)
        colnames(sequences)[1] <- "seq"

        if(fastaExt){
            print("Compressing fasta file")
            gzip(genome)
        }


    } else {
        # Obtain the sequence
        sequences <- data.frame(getSeq(genome, cutSites))
        # Add the cutIDs into the data frame
        sequences$ID <- row.names(sequences)
        colnames(sequences)[1] <- "seq"
    }

    # Subset for cut sites with the desired sequence
    sequences <- sequences[which(sequences$seq == seq),]

    # Turn the sequences into a GRanges object and return it
    sequences <- GRanges(sequences$ID)
    return(sequences)

}
