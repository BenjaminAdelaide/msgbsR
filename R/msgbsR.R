#' msgbsR
#'
#' @name msgbsR
#' @docType package
NULL

#' Read counts of potential MspI cut sites from a MS-GBS experiment of prostates from rats
#'
#' A RangedSummarizedExperiment containing read counts generated from a MS-GBS experiment using the restriction enzyme MspI, focusing on chromosome 20 of Rat.
#'
#' \itemize{
#'		\item ratdata A RangedSummarizedExperiment with 16047 potential MspI cut sites on chromosome 20 in Rat and six samples (3 Control and 3 Experimental).
#' }
#' @docType data
#' @keywords datasets
#' @name ratdata
#' @usage data(ratdata)
#' @format RangedSummarizedExperiment
#' @details This dataset contains six prostate samples from rats: 3 control and 3 experimental high fat diet.
#' @return RangedSummarizedExperiment
NULL

#' Read counts of correct MspI cut sites from a MS-GBS experiment of prostates from rats
#'
#' A RangedSummarizedExperiment containing read counts generated from a MS-GBS experiment using the restriction enzyme MspI, focusing on chromosome 20 of Rat. The sites have been checked for the correct recognition site.
#'
#' \itemize{
#'		\item ratdata2 A RangedSummarizedExperiment containing data for 13983 MspI cut sites on chromosome 20 in Rat and six samples (3 Control and 3 Experimental).
#' }
#' @docType data
#' @keywords datasets
#' @name ratdata2
#' @usage data(ratdata2)
#' @format RangedSummarizedExperiment
#' @details This dataset contains six prostate samples from rats: 3 control and 3 experimental high fat diet. The data can be used for differential methylation analyses.
#' @return RangedSummarizedExperiment
NULL

#' A GRanges object of differentially methylated MspI cut sites on chromosome 20 in Rat from a MS-GBS experiment.
#'
#' The GRanges object was created from a list of differentially methylated cut sites from a MS-GBS experiment between two groups of rats that were fed either a control diet or a high fat diet.
#'
#' \itemize{
#'		\item Positions of MspI cut sites differentially methylated in the prostate on chromosome 20 in Rats.
#' }
#' @docType data
#' @keywords datasets
#' @name cuts
#' @usage data(cuts)
#' @format A GRanges object of length 10.
#' @details The data set contains 10 differentially methylated sites in the prostate between rats fed a control or high fat diet.
#' @return A GRanges object of length 10.
NULL
