#' plotCounts
#'
#' Plot the total number of reads vs total number of cut sites per sample
#'
#' @param countMatrix A data frame containing read counts per sample.
#' @param condition1 A character vector containing a certain condition of the samples.
#' @param condition2 A character vector containing another certain condition of the samples. Default is NULL.
#' @return Produces a plot showing the total number reads vs total number of cut site per sample.
#' @import ggplot2
#' @author Benjamin Mayne
#' @examples
#' data <- system.file("extdata", "datCounts_filtered.Rdata", package = 'msgbsR')
#' load(data)
#' y <- data.frame(c(rep('Control', 3), rep('Fat Diet', 3)))
#' colnames(y) <- 'Group'
#' plotCounts(countMatrix = datCounts, condition1 = y$Group)
#' @export



plotCounts <- function(countMatrix, condition1, condition2 = NULL){

  # Unit tests
  ## Check if the countMatrix is a matrix
  countMatrixClass <- class(countMatrix)[1]
  countMatrixClass <- ifelse(countMatrixClass == 'matrix', yes=TRUE, no=FALSE)
  if(countMatrixClass == 'FALSE'){
    stop('countMatrix must be a matrix')
  }

  # Determine the total number of cuts sites per sample
  # cut sites with > 1 read produced for each sample
  cuts <- data.frame(t(data.frame(lapply(1:ncol(countMatrix),function(x)
    length(countMatrix[,x][which(countMatrix[,x] > 0)])))))

  # Determine the library size (total number of reads)
  libSize <- data.frame(t(data.frame(lapply(1:ncol(countMatrix),
                                            function(x) sum(countMatrix[,x])))))

  if(is.null(condition2)){

    # Make a data frame out of the cuts and libSize
    datPlot <- as.data.frame(cbind(cuts[,1], libSize[,1]))
    datPlot[,3] <- condition1
    colnames(datPlot) <- c('cuts', 'libSize', 'condition1')

    # If else statment depending on whether or not there is a second phenotype characteristic
    qplot(x = libSize, y = cuts, colour = condition1, data = datPlot,
          xlab = 'Total number of Reads per sample',
          ylab = 'Total number of cut sites per sample',
          main = '')

  } else {

  # Make a data frame out of the cuts and libSize
  datPlot <- as.data.frame(cbind(cuts[,1], libSize[,1]))
  datPlot[,3] <- condition1
  datPlot[,4] <- condition2
  colnames(datPlot) <- c('cuts', 'libSize', 'condition1', 'condition2')

  # If else statment depending on whether or not there is a second phenotype characteristic
    qplot(x = libSize, y = cuts, colour = condition1, shape = condition2, data = datPlot,
          xlab = 'Total number of Reads per sample',
          ylab = 'Total number of cut sites per sample',
          main = '')

  }

}
