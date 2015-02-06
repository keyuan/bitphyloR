#' R6 class for inference via MCMC for TSSB.
#'
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @keywords data

TssbMCMC <- R6Class(
  classname = "TssbMCMC",
  inherit = TSSB

  )