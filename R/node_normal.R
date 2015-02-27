#'@include node.R
NULL

ConvertFunctionNameToVariableName <- function(f) {
  res <- unlist(strsplit(gsub("(.)([[:upper:]])", "\\1 \\2", f), " "))
  res[2] <- paste(tolower(substring(res[2], 1, 1)), substring(res[2], 2), sep = "" )
  return(paste(res[-1], collapse = "" ))
}

GetHyperParamsTemplete <- function(f) {
  v <- ConvertFunctionNameToVariableName(f)
  function() {
    print(v)
    if (is.null(self$GetParent())) {
      return(eval(paste("private$", v, sep ="")))
    } else {
      return(do.call(paste("self$GetParent()$", f, sep = ""  ), list()))
    }
  }
}

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
        self$sigma <- riwish(v = priorSigmaDof, S = priorSigmaScale)
        self$params <- rmvnorm(n = 1,
                               mean = initMean,
                               sigma = private$drift)
      } else {
        self$sigma <- riwish(v = priorSigmaDof, S = priorSigmaScale)
        self$params <- rmvnorm(n = 1,
                               mean = parent$params,
                               sigma = self$GetDrift())
      }
    },

    GetDrift = GetHyperParamsTemplete("GetDrift"),

    GetDrift2 = function() {private$drift},

    GetLogProb = function(x) {
      sum(dmvnorm(x, self$params, self$sigma, log = TRUE))
    },

    GetNodeLogProb = function() {
      self$GetLogProb(self$GetData())
    },

    ResampleParams = function() {
      data <- self$GetData()
      drift <- GetDrift


    }
  ),
  private = list(
    drift = NA,
    priorDriftScale = NA,
    priorDriftDof = NA,
    priorSigmaScale = NA,
    priorSigmaDof = NA
  )
)