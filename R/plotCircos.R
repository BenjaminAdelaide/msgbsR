#' plotCircos
#'
#' Plot a circos representing the cut site locations
#'
#' @param cutSites A GRanges object containing the locations of the cut sites to be plotted.
#' @param seqlengths An integer with the lengths of the chromosomes.
#' @param cutSite.colour The colour of the cut sites.
#' @param seqlengths.colour The colour of the chromosomes
#' @usage plotCircos(cutSites, seqlengths, cutSite.colour, seqlengths.colour)
#' @return A circos plot showing the locations of the cut sites.
#' @importFrom  S4Vectors Rle
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @importFrom IRanges IRanges
#' @importFrom ggbio autoplot circle ggbio
#' @import ggplot2
#' @author Benjamin Mayne
#' @examples
#' # load example cut site positions
#' data(cuts)
#' # Obtain the length of chromosome 20 in rn6
#' library(BSgenome.Rnorvegicus.UCSC.rn6)
#' chr20 <- seqlengths(BSgenome.Rnorvegicus.UCSC.rn6)["chr20"]
#' plotCircos(cutSites = cuts, seqlengths = chr20,
#'            cutSite.colour = "red", seqlengths.colour = "blue")
#' @export

plotCircos <- function(cutSites, seqlengths, cutSite.colour, seqlengths.colour){

  # Make a matrix out of the genome object to turn into a GRanges object
  seqlengths <- as.matrix(seqlengths)

  # Add a column into the seqlengths matrix of 1s and then reshape it for a GRanges object
  seqlengths <- data.frame(as.matrix(cbind(seqlengths, rep(1, nrow(seqlengths)))))
  seqlengths[,3] <- row.names(seqlengths)
  seqlengths <- as.matrix(data.frame(seqlengths)[ ,c(3,2,1)])
  seqlengths.gr <- GRanges(seqnames = Rle(values = c(seqlengths[,1]), lengths = as.numeric(rep("1", nrow(seqlengths)))),
                       ranges = IRanges(start = as.numeric(seqlengths[,2]), end = as.numeric(seqlengths[,3])))

  seqlevels(seqlengths.gr) <- as.character(seqlengths[,1])
  seqlengths(seqlengths.gr) <- as.integer(seqlengths[,3])

  seqlevels(cutSites) = seqlevels(seqlengths.gr)
  seqlengths(cutSites)=as.integer(seqlengths[,3])

  # Plot using specified colours
  p <- ggbio() + circle(cutSites, geom = "rect", color = cutSite.colour) +
    circle(seqlengths.gr, geom = "ideo", fill = seqlengths.colour) +
    circle(seqlengths.gr, geom = "scale", size = 2) +
    circle(seqlengths.gr, geom = "text", aes(label = seqnames), vjust = 0, size = 3)
  p

}
