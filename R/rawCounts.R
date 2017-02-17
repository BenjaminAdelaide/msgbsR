#' rawCounts
#'
#' Imports the raw read counts from sorted and indexed bam file(s)
#'
#' @param bamFilepath The path to the location of the bam file(s).
#' @param threads The total number of usable threads to be used. Default is 1.
#' @usage rawCounts(bamFilepath, threads = 1)
#' @return Produces a RangedSummarizedExperiment. Columns are samples and the rows are cut sites. The cut site IDs are in the format chr:position-position:strand.
#' @importFrom easyRNASeq validate BamFileList
#' @import Rsamtools
#' @import GenomicFeatures
#' @import GenomicRanges
#' @import GenomicAlignments
#' @import plyr
#' @import parallel
#' @import SummarizedExperiment
#' @author Benjamin Mayne, Sam Buckberry
#' @examples
#'my_path <- system.file("extdata", package = "msgbsR")
#'my_data <- rawCounts(bamFilepath = my_path)
#' @export


# Function to make the raw count matrix
rawCounts <- function(bamFilepath, threads = 1){

    # function to get counts for start of each read
    readCounts <- function(bamFiles){
        # Set the features to extract for each BAM record
        what <- c("rname", "strand", "pos", "qwidth")
        # Set the scan parameters
        scanParam <- ScanBamParam(flag=scanBamFlag(isDuplicate=NA,
                                                   isUnmappedQuery=FALSE),
                                  what=what)

        gal <- readGAlignments(file = bamFiles, param = scanParam)

        # Get the read starts
        readStart <- ifelse(strand(gal) == "-", end(gal), start(gal))

        # convert to a data frame
        df <- data.frame(gal)

        # add readStart
        df$start <- readStart

        # Count the number of reads starting at each locus
        out <- ddply(df, c("rname", "strand", "start"), .fun=nrow)
        colnames(out)[4] <- "count"

        # Create a unique locus identifier
        out$id <- paste(out$rname, out$start, sep=":")
        out$id <- paste(out$id, out$start, sep="-")
        out$id <- paste(out$id, out$strand, sep=":")

        # Add the sample name
        out$sample <- basename(bamFiles)

        # output the dataframe
        out
    }

    #============================================================================================
    # Firstly check if the threads input vaue is numeric
    if(!is(threads, "numeric")){
        stop("The threads must be a numeric value")
    }

    # Get a list of all the bam files
    bamFiles <- list.files(path=bamFilepath,
                           full.names=TRUE, pattern=".bam$")

    # Get a list of all the indexed files (*.bai)
    baiFiles <- list.files(path=bamFilepath,
                           full.names=TRUE, pattern=".bai$")

    # Make a stop function if a bam file is missing an indexed file
    BAMS <- basename(bamFiles)
    BAIS <- gsub(".bai", "", basename(baiFiles))

    NoBAIS <- BAMS[which(BAMS %in% BAIS == 'FALSE')]

    if(length(NoBAIS) != 0){
        message <- paste(c("The following sample(s) are missing indexed files", NoBAIS), sep='\t', collapse = ', ')
        message <- gsub("files,", "files:", message)
        stop(message)
    }

    # Check that the BAM files are missing a EOF header
    validate(BamFileList(bamFiles,index=bamFiles))

    # Run readCounts on each bam with multiple threads
    dat <- mclapply(X = bamFiles, FUN=readCounts, mc.cores=threads)

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

        #First, get the location of the matching ids for the sample in the matrix
        index <- match(x=sampleDat$id, table=row.names(countMatrix))

        counts <- sampleDat$count

        countMatrix[index, i] <- counts
    }

    # Remove the .bam extension from the colNames
    colnames(countMatrix) <- gsub(pattern=".bam", replacement="", colnames(countMatrix))

    # Convert to an integer matrix
    countMatrix <- apply(countMatrix, c(1, 2), function(x) {(as.integer(x))})

    # Output the data as a RangedSummarizedExperiment
    se <- SummarizedExperiment(assays=list(counts=countMatrix),
                         rowRanges=GRanges(rownames(countMatrix)), colData=colnames(countMatrix))

    return(se)
}








