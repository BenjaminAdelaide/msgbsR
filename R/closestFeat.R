#' closestFeat
#'
#' Determines the closest feature of cut site locations
#'
#' @param gff The path to the location of a gff file.
#' @param cutSites A GRanges object containing the location of the cut sites.
#' @return The closest biological feature annotated in a gff3 file to each cut site.
#' @author Benjamin Mayne
#' @import GenomicRanges
#' @importFrom utils read.table
#' @examples
#' data(cuts)
#' mygff <- system.file("extdata", "chr20.gff3", package = "msgbsR")
#' features <- closestFeat(gff = mygff, cutSites = cuts)
#' @export

closestFeat <- function(gff, cutSites){

  ## Check if the gff is a gff file
  extenstion <- function (x)
  {
    pos <- regexpr("\\.([[:alnum:]]+)$", x)
    ifelse(pos > -1L, substring(x, pos + 1L), "")
  }
  gffExt <- ifelse(extenstion(gff) == "gff" | extenstion(gff) == "gff3" , yes=TRUE, no=FALSE)
  if(!isTRUE(gffExt)){
    stop("Input is not a gff or gff3 file")
  }
  gffRead <- function(gffFile, nrows = -1) {
    gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
                     header=FALSE, comment.char="#", nrows = nrows,
                     colClasses=c("character", "character", "character", "integer",
                                  "integer",
                                  "character", "character", "character", "character"))
    colnames(gff) = c("seqname", "source", "feature", "start", "end",
                      "score", "strand", "frame", "attributes")

    return(gff)
  }

  # Message
  print("Reading in the gff file")
  # Read in the gff file
  gff <- gffRead(gffFile = gff)

  # Turn into a GRanges object
  gff2 <- GRanges(seqnames = Rle(values = c(gff[,1]), lengths = as.numeric(rep("1", nrow(gff)))),
          ranges = IRanges(start = as.numeric(gff[,4]), end = as.numeric(gff[,5])))
  # Find the nearest feature
  nearestFeature <- gff[nearest(cutSites, gff2),]
  row.names(nearestFeature) <- NULL
  return(nearestFeature)
}

