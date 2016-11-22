#' plotChr
#'
#' Plot a karogram representing the cut site locations
#'
#' @param cutSites A matrix consisting of three columns, where the first column is the chromosome ID, the second and third column being the start and end position of the cut site respectively.
#' @param genome A matrix with the lengths of the chromosomes. The first column is the chromosome ID and the second is the length of each chromosome.
#' @usage plotChr(cutSites, genome)
#' @return An ideogram showing the locations of the cut sites.
#' @importFrom  S4Vectors Rle
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @importFrom IRanges IRanges
#' @importFrom ggbio autoplot circle ggbio layout_karyogram autoplot
#' @import ggplot2
#' @author Benjamin Mayne
#' @examples
#' chr20 <- matrix(c('chr20', '56205956'), nrow = 1, ncol = 2)
#' data <- system.file("extdata", "exampleSites.Rdata", package = 'msgbsR')
#' load(data)
#' plotChr(cutSites = z[,c(1,3,3)], genome = chr20)
#' @export


plotChr <- function(cutSites, genome){

  # Add another column into the genome matrix of 1s
  genome <- as.matrix(cbind(genome, rep(1, nrow(genome))))
  genome <- as.matrix(data.frame(genome)[ ,c(1,3,2)])

  # Turn the data frames into GRanges objects
  cutSites.gr <- GRanges(seqnames = Rle(values = c(cutSites[,1]), lengths = as.numeric(rep('1', nrow(cutSites)))),
                         ranges = IRanges(start = as.numeric(cutSites[,2]), end = as.numeric(cutSites[,3])))
  genome.gr <- GRanges(seqnames = Rle(values = c(genome[,1]), lengths = as.numeric(rep('1', nrow(genome)))),
                       ranges = IRanges(start = as.numeric(genome[,2]), end = as.numeric(genome[,3])))

  #seqnames(cutSites.gr) <- seqnames(genome.gr)
  seqlevels(genome.gr) <- as.character(genome[,1])
  seqlengths(genome.gr) <- as.integer(genome[,3])

  seqlevels(cutSites.gr, force = TRUE) = seqlevels(genome.gr)
  seqlengths(cutSites.gr)=as.integer(genome[,3])

  # Plot the karyogram
  autoplot(seqinfo(genome.gr))
  autoplot(seqinfo(cutSites.gr)) + layout_karyogram(cutSites.gr)

}
