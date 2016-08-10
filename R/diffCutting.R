#' diffCutting
#' 
#' Determines differential cutting from a matrix of read counts
#' 
#' @param countMatrix A data frame containing read counts. Rows are cut sites and columns are samples.
#' @param pd A data frame containing any meta data of the samples within the countMatrix.
#' @param cateogory The column name in pd that is to be tested for differential cutting.
#' @param condition1 The reference group within the cateogory of pd.
#' @param condition2 The experimental group within the cateogory of pd.
#' @param block The column name of pd if differential cutting is to be tested with a blocking factor. Default is NULL.
#' @param cpmThreshold Counts per million threshold of read counts to be filtered out of the analysis.
#' @param thresholdSamples Minimum number of samples to contain the counts per million threshold.
#' @return A data frame containing which cut sites are differenitally cut.
#' @examples
#' \dontrun{
#' 
#' # Make a read counts data frame 
#' set.seed(1)
#' x <- data.frame(matrix(sample(0:100, size = 10000*10, replace = TRUE), nrow = 10000, ncol = 10))
#' # Set the groups and blocking factor
#' y <- data.frame(c(rep('A', 5), rep('B', 5)))
#' colnames(y) <- 'group'
#' y$Block <- rep(c(rep('C',1), rep('D', 1)), 5)
#' 
#' # Determine differential cut sites
#' z <- diffCutting(countMatrix = x, pd = y, cateogory = 'group', condition1 = 'A',
#'               condition2 = 'B', block = 'Block', cpmThreshold = 1, thresholdSamples =3)
#' 
#' }
#' @author Benjamin Mayne
#' @export


# Function for determining differential cutting

diffCutting <- function(countMatrix, pd, cateogory, condition1, condition2, block = NULL,
                        cpmThreshold, thresholdSamples){
  
  # Seperate the data into the two conditions
  counts <- data.frame(countMatrix[ ,pd[,cateogory] == condition1 | pd[,cateogory] == condition2])
  group <- pd[,cateogory][pd[,cateogory] == condition1 | pd[,cateogory] == condition2]
  pd <- pd[pd[,cateogory] == condition1 | pd[,cateogory] == condition2,]
  
  # Make the DGE object
  dge <- DGEList(counts=counts, 
                 group=group, 
                 genes = row.names(counts))
  
  # Remove lowly cut sites as defined by user
  dge_expressed <- dge[rowSums(cpm(dge) > cpmThreshold) >= thresholdSamples, ]
  
  # Re-compute the library sizes after filtering
  dge_expressed$samples$lib.size <- colSums(dge_expressed$counts)
  
  # Data Normalisation
  # Compute effective library sizes using the default trimmed mean of M-values (TMM)   
  dge_expressed <- calcNormFactors(dge_expressed, method = "TMM")
  
  if(is.null(block)){

    group = factor(group)
    group <- relevel(group, ref = condition1)
    design <- model.matrix(~group)
    
  } else {
    
    # Make the design matrix and use the blocking
    block = factor(pd[,block])
    group = factor(group)
    group <- relevel(group, ref = condition1)
    design <- model.matrix(~block+group)
    
  }
  
  # Estimate the dispersion parameters and fit the model
  dispersions <- estimateDisp(dge_expressed, design)
  
  # fit the genewise negative binomial GLMs for design
  fit <- glmFit(dispersions, design)
  
  # Statistical Analyses for differential methylaiton 
  # Coefficient contrast
  if(is.null(block)){
    lrt <- glmLRT(fit, coef = colnames(design)[2])
  } else {
    lrt <- glmLRT(fit, coef = colnames(design)[3])
  }  
  
  # Use Benjamini-Hochberg (BH) method for p-values
  lrt_top <- topTags(lrt, n = nrow(dge_expressed$counts), adjust.method = "BH", sort.by = "PValue")
  
  return(lrt_top$table)
}






