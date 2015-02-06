#' R6 class for inference via stochastic variational Bayes for TSSB.
#'
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @keywords data

TssbSVB <- R6Class(
  classname = "TssbSVB",
  inherit = TSSB

)