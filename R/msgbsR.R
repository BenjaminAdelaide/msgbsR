#' msgbsR
#'
#' @name msgbsR
#' @docType package
NULL

#' Table of read counts
#'
#' A data set containing the following data:
#'
#' \itemize{
#'		\item datCounts a matrix with 16977 rows (Potential MspI cut sites) and six columns (samples)
#' }
#' @docType data
#' @keywords datasets
#' @name datCounts
#' @usage data(datCounts)
#' @format A 16977x6 matrix
#' @return A 16977x6 matrix
NULL

#' Table of filtered read counts from a MSGBS experiment
#'
#' A data set containing the following data:
#'
#' \itemize{
#'		\item datCounts_filtered a matrix with 7178 rows (Correct MspI cut sites) and six columns (samples)
#' }
#' @docType data
#' @keywords datasets
#' @name datCounts_filtered
#' @usage data(datCounts_filtered)
#' @format A 7178x6 matrix
#' @return A 7178x6 matrix
NULL

#' A integer with the length of chr20 in Rat
#'
#' A data set containing the following data:
#'
#' \itemize{
#'		\item ratChr an integer of the value 56205956
#' }
#' @docType data
#' @keywords datasets
#' @name ratChr
#' @usage data(ratChr)
#' @format A integer of length 1
#' @return A integer of length 1
NULL


#' A GRanges object of MspI cut sites.
#'
#' A data set containing the following data:
#'
#' \itemize{
#'		\item cuts a GRanges object of MspI cut sites
#' }
#' @docType data
#' @keywords datasets
#' @name cuts
#' @usage data(cuts)
#' @format A GRanges object
#' @return A GRanges object
NULL
