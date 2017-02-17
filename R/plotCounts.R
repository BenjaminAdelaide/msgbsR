#' plotCounts
#'
#' Plots the total number of reads vs total number of cut sites per sample
#'
#' @param se A RangedSummarizedExperiment containing meta data of the samples.
#' @param cateogory The heading name in the sample data to distinguish groups.
#' @return Produces a plot showing the total number reads vs total number of cut sites per sample.
#' @usage plotCounts(se, cateogory)
#' @import ggplot2
#' @author Benjamin Mayne
#' @examples
#' data(ratdata2)
#' plotCounts(se = ratdata2, cateogory = "Group")
#' @export

plotCounts <- function(se, cateogory){

    # Unit tests
    ## Check if input is a RangedSummarizedExperiment
    if(!is(se, "RangedSummarizedExperiment")){
    stop("se must be a RangedSummarizedExperiment")
    }

    ## Check if cateogory is a character
    if(!is(cateogory, "character")){
        stop("cateogory must be a character")
    }

    # Determine the total number of cuts sites per sample
    # cut sites with > 1 read produced for each sample
    cuts <- data.frame(t(data.frame(lapply(1:ncol(assay(se)),function(x)
    length(assay(se)[,x][which(assay(se)[,x] > 0)])))))

    # Determine the library size (total number of reads)
    libSize <- data.frame(t(data.frame(lapply(1:ncol(assay(se)),
                                            function(x) sum(assay(se)[,x])))))

    # Make a data frame out of the cuts and libSize
    datPlot <- as.data.frame(cbind(cuts[,1], libSize[,1]))
    datPlot[,3] <- colData(se)[,cateogory]
    colnames(datPlot) <- c("cuts", "libSize", "cateogory")

    # Plot the total number of counts vs total number of cuts
    qplot(x = libSize, y = cuts, colour = cateogory, data = datPlot,
          xlab = "Total number of Reads per sample",
          ylab = "Total number of cut sites per sample",
          main = "")

}
