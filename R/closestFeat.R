#' closestFeat
#'
#' Determines the closest feature of cut site locations
#'
#' @param gff The path to the location of a gff file.
#' @param positions A matirx where the first column is the chromosome ID, the second column is the strand (1 or 2) and the third column is the position.
#' @param strand.specific If TRUE the function will return the closest feature on the same strand, otherwise it'll return the closest feature on either strand. Default is TRUE.
#' @author Benjamin Mayne
#' @export

closestFeat <- function(gff, positions, strand.specific = TRUE){

  # Read in the gff file
  gff <- gffRead(gffFile = gff)

  # Message
  cat("Finding closest feature...", "\n")

  # Start of looping through each position
  theclosestFeatures <- lapply(1:nrow(positions), function(x) findClosest(gff = gff,
              positions = positions, strand.specific = strand.specific, x = x))
  # Unlist the object
  theclosestFeatures <- ldply(theclosestFeatures)

  return(theclosestFeatures)
}

gffRead <- function(gffFile, nrows = -1) {
  cat("Reading ", gffFile, ": ", sep="")
  gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
                   header=FALSE, comment.char="#", nrows = nrows,
                   colClasses=c("character", "character", "character", "integer",
                                "integer",
                                "character", "character", "character", "character"))
  colnames(gff) = c("seqname", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attributes")
  return(gff)
}

findClosest <- function(gff, positions, strand.specific = TRUE, x){

  # Seperate the gff file into the chromosome of interest
  gffChr <- gff[which(gff$seqname == positions[x,1]), ]

  # Determine the closest feature in the gffChr file
  if(strand.specific == TRUE){

    # Match up the strands (1='+' 2='-')
    strand <- as.character(ifelse(test = positions[x,2] == '1', yes = '+', no = '-'))
    gffChr <- gffChr[which(gffChr$strand == strand), ]

  }

  # Determine the distance from each start position
  dist2start <- data.frame(as.numeric(rep(positions[x,3], nrow(gffChr))) - gffChr$start)
  dist2start[,2] <- abs(as.numeric(dist2start[,1]))
  rowNo <- which(dist2start[,2] == min(dist2start[,2])[1])[1]

  theFeature <- data.frame(cbind(dist2start[rowNo, 1], gffChr[rowNo,c(1,3:9)]))
  colnames(theFeature)[1] <- 'Distance'
  return(theFeature)

}














