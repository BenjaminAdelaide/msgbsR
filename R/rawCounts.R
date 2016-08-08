#' rawCounts
#'
#' Imports the raw read counts from sorted and indexed bam file(s)
#'
#' @param bamFilePath The path to the location of bam file(s).
#' @param threads The total number of usable threads to be used. Default is 1.
#' @usage datCounts <- rawCounts(bamFilepath = 'path-to-bam-files', threads = 1)
#' @return Produces a data frame where the rows correspond to cut sites and columns are the samples. Cut sites are given a unique ID (chromosome:strand:position).
#' @author Benjamin Mayne, Sam Buckberry
#' @export

# Function to make the raw count matrix
rawCounts <- function(bamFilepath, threads = 1){

  # function to get counts for start of each read
  readCounts <- function(bamFilePath){
    # Set the features to extract for each BAM record
    what <- c("rname", "strand", "pos", "qwidth")
    # Set the scan parameters
    scanParam <- ScanBamParam(flag=scanBamFlag(isDuplicate=NA,
                                               isUnmappedQuery=FALSE),
                              what=what)

    bam <- scanBam(bamFilePath, param=scanParam)

    # Get the bam data in a list as per the Rsamtools documentation
    lst <- lapply(names(bam[[1]]), function(elt) {
      do.call(c, unname(lapply(bam, "[[", elt)))
    })
    names(lst) <- names(bam[[1]])

    df <- do.call("DataFrame", lst)

    # Get the position of the read for + and - strand. Remember the mapping
    # to the minus strand has the read in the opposite orientation
    # Minus 1 for the zero based coordinates
    df$readStart <- as.integer(ifelse(test=df$strand == 2,
                                      yes=df$pos + df$qwidth - 1, no=df$pos))

    # Set the chromosome names
    rname_new <- bam[[1]][1]
    rname_new <- as.character(rname_new$rname)
    df$rname <- rname_new

    # remove the irrelevant columns
    df <- df[ ,c(1, 2, 5)]
    df <- data.frame(df)
    head(df)

    # Count the number of reads starting at each locus
    out <- ddply(df, .(rname, strand, readStart), .fun=nrow)
    colnames(out)[4] <- "count"

    # Create a unique locus identifier
    out$id <- paste(out$rname, out$strand, out$readStart, sep=":")

    # Add the sample name
    out$sample <- basename(bamFilePath)

    # output the dataframe
    out
  }


  # Get a list of all the bam files
  bamFiles <- list.files(path=bamFilepath,
                         full.names=TRUE, pattern=".bam$")

  # Use the function readCounts to get the start position of each read
  dat <- readCounts(bamFilePath=bamFiles[1])

  ###### Insert a stop function here if any are duplicated ##########
  anyDuplicated(dat$id)

  # Run readCounts on each bam with multiple threads
  dat <- mclapply(bamFiles, FUN=readCounts, mc.cores=threads)

  # Turn the list into a data frame
  dat <- ldply(dat)

  # Get all the unique cut sites
  ids <- unique(dat$id)

  # Setup the raw count matrix
  countMatrix <- matrix(data=0, nrow=length(ids), ncol=length(bamFiles))
  row.names(countMatrix) <- ids
  colnames(countMatrix) <- basename(bamFiles)

  # Function to match locus ids in each sample and get counts
  bams <- basename(bamFiles)

  for(i in 1:length(bams)){
    # Subset the sample of interest
    sampleDat <- dat[dat$sample == bams[i], ]

    # Test for unique locus ids and return error if false
    stopifnot(anyDuplicated(sampleDat$id) == 0)

    #First, get the location of the matching ids for the sample in the matrix
    index <- match(x=sampleDat$id, table=row.names(countMatrix))

    counts <- sampleDat$count

    countMatrix[index, i] <- counts
  }

  # Remove the .bam extension from the colNames
  colnames(countMatrix) <- gsub(pattern=".bam", replacement="", colnames(countMatrix))

  return(countMatrix)
}
