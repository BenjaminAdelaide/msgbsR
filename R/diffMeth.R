#' diffMeth
#'
#' Determines differential methylated sites from a RangedSummarizedExperiment
#'
#' @param se A RangedSummarizedExperiment containing meta data of the samples.
#' @param cateogory The heading name in the sample data to be tested for differential methylation.
#' @param condition1 The reference group within the cateogory.
#' @param condition2 The experimental group within the cateogory.
#' @param block The heading name in the sample data if differential methylation is to be tested with a blocking factor. Default is NULL.
#' @param cpmThreshold Counts per million threshold of read counts to be filtered out of the analysis.
#' @param thresholdSamples Minimum number of samples to contain the counts per million threshold.
#' @usage diffMeth(se, cateogory, condition1, condition2,
#'                  block = NULL, cpmThreshold, thresholdSamples)
#' @return A data frame containing which cut sites that are differenitally methylated.
#' @import edgeR
#' @importFrom stats model.matrix relevel
#' @author Benjamin Mayne
#' @examples
#' # Load data
#' data(ratdata2)
#' top <- diffMeth(se = ratdata2, cateogory = "Group",
#'        condition1 = "Control", condition2 = "Experimental",
#'        cpmThreshold = 1, thresholdSamples = 1)
#' @export

# Function for determining differential cutting/methylation

diffMeth <- function(se, cateogory, condition1, condition2, block = NULL,
                        cpmThreshold, thresholdSamples){

    # Unit tests
    ## Check if se is a RangedSummarizedExperiment
    if(!is(se, "RangedSummarizedExperiment")){
    stop("se must be a RangedSummarizedExperiment")
    }

    ## Check if cateogory is a character
    if(!is(cateogory, "character")){
        stop("cateogory must be a character")
    }

    ## Check if condition1 is a character
    if(!is(condition1, "character")){
    stop("condition1 must be a character")
    }

    ## Check if condition2 is a character
    if(!is(condition2, "character")){
    stop("condition2 must be a character")
    }

    ## Check if cpmThreshold is a numeric
    if(!is(cpmThreshold, "numeric")){
    stop("cpmThreshold must be a numeric")
    }

    ## Check if thresholdSamples is a numeric
    if(!is(thresholdSamples, "numeric")){
    stop("thresholdSamples must be a numeric")
    }

    pd <- colData(se)

    # Seperate the data into the two conditions
    counts <- data.frame(assay(se)[ ,pd[,cateogory] == condition1 | pd[,cateogory] == condition2])
    group <- pd[,cateogory][pd[,cateogory] == condition1 | pd[,cateogory] == condition2]
    pd <- pd[pd[,cateogory] == condition1 | pd[,cateogory] == condition2,]

    # Make the DGE object
    dge <- DGEList(counts=counts, group=group,
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
    colnames(lrt_top$table)[1] <- "site"
    row.names(lrt_top$table) <- NULL
    return(lrt_top$table)

}

