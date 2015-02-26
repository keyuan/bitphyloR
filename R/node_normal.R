#'@include node.R
NULL

#' R6 class for Normal node.
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @keywords Normal node class
#' @field dataIds: data IDs
#' @field tssb: A TSSB object
#' @field children: descentant node objects
#' @field parent: parent node object
#' @method initialize
#' @method GetChildren
Normal <- R6::R6Class(
  classname = "Normal",
  inherit = Node,
  public = list(
    initialize = function(parent = NULL, tssb = NULL) {
      super$new(parent = parent, tssb = tssb)
      if (is.null(parent)) {
        private$drift <- drift
        private$priorDriftScale <- priorDriftScale
        private$priorDriftDof <-priorDriftDof
        private$priorSigmaScale <- priorSigmaScale
        private$priorSigmaDof <- priorSigmaDof
        self$sigma <- MCMCpack::riwish(v = priorSigmaDof, S = priorSigmaScale)
        self$params <- mvtnorm::rmvnorm(n = 1,
                                        mean = initMean,
                                        sigma = private$drift)
      } else {
        self$sigma <- MCMCpack::riwish(v = priorSigmaDof, S = priorSigmaScale)
        self$params <- mvtnorm::rmvnorm(n = 1,
                                        mean = parent$params,
                                        sigma = self$GetDrift())
      }

    }


    ),
  private = list(
    drift = NA,
    priorDriftScale = NA,
    priorDriftDof = NA
    )
)