#' plotCircos
#'
#' Plot a circos representing the cut site locations
#'
#' @param cutSites A matrix consisting of three columns, where the first column is the chromosome ID, the second and third column being the start and end position of the cut site respectively.
#' @param genome An integer with the lengths of the chromosomes.
#' @param cutSite.colour The colour of the cut sites.
#' @param genome.colour The colour of the chromosomes
#' @usage plotCircos(cutSites, genome, cutSite.colour, genome.colour)
#' @return A circos plot showing the locations of the cut sites.
#' @importFrom  S4Vectors Rle
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @importFrom IRanges IRanges
#' @importFrom ggbio autoplot circle ggbio
#' @import ggplot2
#' @author Benjamin Mayne
#' @examples
#' data(cuts)
#' data(ratChr)
#' plotCircos(cutSites = cuts, genome = ratChr,
#'            cutSite.colour = 'red', genome.colour = 'blue')
#' @export


plotCircos <- function(cutSites, genome, cutSite.colour, genome.colour){

  # Make a matrix out of the genome object to turn into a GRanges object
  genome <- as.matrix(genome)

  # Add a column into the genome matrix of 1s and then reshape it for a GRanges object
  genome <- data.frame(as.matrix(cbind(genome, rep(1, nrow(genome)))))
  genome[,3] <- row.names(genome)
  genome <- as.matrix(data.frame(genome)[ ,c(3,2,1)])
  genome.gr <- GRanges(seqnames = Rle(values = c(genome[,1]), lengths = as.numeric(rep('1', nrow(genome)))),
                       ranges = IRanges(start = as.numeric(genome[,2]), end = as.numeric(genome[,3])))

  seqlevels(genome.gr) <- as.character(genome[,1])
  seqlengths(genome.gr) <- as.integer(genome[,3])

  seqlevels(cutSites, force = TRUE) = seqlevels(genome.gr)
  seqlengths(cutSites)=as.integer(genome[,3])

  # Plot using specified colours
  p <- ggbio() + circle(cutSites, geom = "rect", color = cutSite.colour) +
    circle(genome.gr, geom = "ideo", fill = genome.colour) +
    circle(genome.gr, geom = "scale", size = 2) +
    circle(genome.gr, geom = "text", aes(label = seqnames), vjust = 0, size = 3)
  p

}
