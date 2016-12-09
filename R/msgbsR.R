#' msgbsR
#'
#' @name msgbsR
#' @docType package
NULL

#' A table of read counts from a MS-GBS experiment
#'
#' An integer matrix containing read counts generated from a MS-GBS experiment using the restriction enzyme MspI, focusing on chromosome 20 of Rat.
#'
#' \itemize{
#'		\item datCounts a matrix with 16977 rows (Potential MspI cut sites on chromosome 20) and six columns (samples).
#' }
#' @docType data
#' @keywords datasets
#' @name datCounts
#' @usage data(datCounts)
#' @format A 16977x6 integer matrix
#' @return A 16977x6 integer matrix
NULL

#' A table of read counts from a MS-GBS experiment after filtering incorrect cut sites
#'
#' An integer matrix of read counts with correct locations of MspI cuts of chromosome 20 from Rat.
#'
#' \itemize{
#'		\item datCounts_filtered a matrix with 7178 rows (Correct MspI cut sites) and six columns (samples)
#' }
#' @docType data
#' @keywords datasets
#' @name datCounts_filtered
#' @usage data(datCounts_filtered)
#' @format A 7178x6 integer matrix
#' @return A 7178x6 integer matrix
NULL

#' The length of chromsome 20 in Rats (rn6 UCSC)
#'
#' An integer with the length of base pairs of chr20 in Rat (rn6 UCSC).
#'
#' \itemize{
#'		\item ratChr an integer of the value 56205956
#' }
#' @docType data
#' @keywords datasets
#' @name ratChr
#' @usage data(ratChr)
#' @format An integer of length 1
#' @return An integer of length 1
NULL

#' A GRanges object of differentially methylated MspI cut sites from a MS-GBS experiment.
#'
#' The GRanges object was created from a list of differentially methylated cut sites from a MS-GBS experiment between two groups of rats that were fed different diets.
#'
#' \itemize{
#'		\item Positions of MspI cut sites differentially methylated on chromosome 20 in Rats.
#' }
#' @docType data
#' @keywords datasets
#' @name cuts
#' @usage data(cuts)
#' @format A GRanges object of length 363.
#' @return A GRanges object of length 363.
NULL
