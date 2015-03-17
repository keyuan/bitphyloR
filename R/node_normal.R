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
    sigma = NULL,
    params = NULL,

    initialize = function(parent = NULL, tssb = NULL,
                          drift = diag(1), priorDriftScale = diag(1), priorDriftDof = 1,
                          priorSigmaScale = diag(1), priorSigmaDof = 1, initMean = 0) {
      super$initialize(parent = parent, tssb = tssb)
      if (is.null(parent)) {
        private$drift <- drift
        private$priorDriftScale <- priorDriftScale
        private$priorDriftDof <-priorDriftDof
        private$priorSigmaScale <- priorSigmaScale
        private$priorSigmaDof <- priorSigmaDof
        private$initMean <- initMean
        self$sigma <- riwish(v = priorSigmaDof, S = priorSigmaScale)
        self$params <- rmvnorm(n = 1,
                               mean = initMean,
                               sigma = private$drift)
      } else {
        self$sigma <- riwish(v = priorSigmaDof, S = priorSigmaScale)
        self$params <- rmvnorm(n = 1,
                               mean = parent$params,
                               sigma = parent$drift)
      }
    },

    GetDrift = function() {
      if (is.null(private$parent)) {
        private$drift
      } else {
        private$parent$GetDrift()
      }
    },

    GetLogProb = function(x) {
      sum(dmvnorm(x, self$params, self$sigma, log = TRUE))
    },

    GetNodeLogProb = function() {
      self$GetLogProb(self$GetData())
    },

    ResampleParams = function() {
      data <- self$GetData()
      drift <- self$GetDrift()
      numOfData <- nrow(data)
      numOfChildren <- length(self$GetChildren())

      if (is.null(private$parent)) {
        parentParams <- self$initMean
      } else {
        parentParams <- private$parent$params
      }

    },

    ResampleHyperParams = function() {

    }
  ),
  private = list(
    drift = NA,
    priorDriftScale = NA,
    priorDriftDof = NA,
    priorSigmaScale = NA,
    priorSigmaDof = NA,
    initMean = NA
  )
)