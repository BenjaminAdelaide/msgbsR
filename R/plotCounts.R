#' plotCounts
#'
#' Plot the total number of reads vs total number of cut sites per sample
#'
#' @param countMatrix A data frame containing read counts per sample.
#' @param condition1 A character vector containing a certain condition of the samples.
#' @param condition2 A character vector containing another certain condition of the samples. Default is NULL.
#' @return Produces a plot showing the total number reads vs total number of cut site per sample.
#' @importFrom ggplot2 qplot
#' @author Benjamin Mayne
#' @examples
#' # Generate a random matrix
#' set.seed(1)
#' x <- data.frame(matrix(sample(0:100, size = 10000*10, replace = TRUE), nrow = 10000, ncol = 10))
#' y <- c(rep('A', 5), rep('B', 5))
#' z <- c(rep('C', 3), rep('D', 2), rep('E', 2), rep('F', 3))
#'
#' plotCounts(x,y,z)
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
