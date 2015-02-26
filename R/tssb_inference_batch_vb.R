#' @include tssb.R
NULL

#' R6 class for inference via Batch variational Bayes for TSSB.
#'
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @keywords data
TssbBatchVB <- R6::R6Class(
  classname = "TssbBatchVB",
  inherit = TSSB
)
