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
                          drift = diag(1), priorDriftScale = diag(1),
                          priorDriftDof = nrow(priorDriftScale) + 1,
                          priorSigmaScale = diag(1),
                          priorSigmaDof = nrow(priorSigmaScale) + 1,
                          initMean = array(0, dim=c(1, nrow(priorSigmaScale)))) {
      super$initialize(parent = parent, tssb = tssb)
      if (is.null(parent)) {
        private$drift <- riwish(v = priorDriftDof, S = priorDriftScale)
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
        self$sigma <- riwish(v = self$GetPriorSigmaDof(),
                             S = self$GetPriorSigmaScale())
        self$params <- rmvnorm(n = 1,
                               mean = parent$params,
                               sigma = self$GetDrift())
      }
    },

    GetPriorDriftDof = function() {
      if (is.null(private$parent)) {
        private$priorDriftDof
      } else {
        private$parent$GetPriorDriftDof()
      }
    },

    GetPriorDriftScale = function() {
      if (is.null(private$parent)) {
        private$priorDriftScale
      } else {
        private$parent$GetPriorDriftScale()
      }
    },

    GetDrift = function() {
      if (is.null(private$parent)) {
        private$drift
      } else {
        private$parent$GetDrift()
      }
    },

    GetPriorSigmaDof = function() {
      if (is.null(private$parent)) {
        private$priorSigmaDof
      } else {
        private$parent$GetPriorSigmaDof()
      }
    },

    GetPriorSigmaScale = function() {
      if (is.null(private$parent)) {
        private$priorSigmaScale
      } else {
        private$parent$GetPriorSigmaScale()
      }
    },

    GetLogProb = function(x) {
      sum(dmvnorm(x, self$params, self$sigma, log = TRUE))
    },

    GetNodeLogProb = function() {
      self$GetLogProb(self$GetData())
    },

    ResampleParams = function() {
      nodeData <- self$GetData()
      drift <- self$GetDrift()
      numOfData <- self$GetNumOfLocalData()
      priorSigmaScale <- self$GetPriorSigmaScale()
      priorSigmaDof <- self$GetPriorSigmaDof()
      numOfChildren <- length(private$children)

      if (numOfData == 0) {
        dataMean = 0
      } else if (numOfData == 1){
        dataMean = nodeData
      } else {
        dataMean = colMeans(nodeData)
      }

      if (is.null(private$parent)) {
        parentParams <- private$initMean
      } else {
        parentParams <- private$parent$params
      }

      if (numOfChildren ==0 ) {
        childParamsMean = 0
      } else {
        childParams <- Reduce(
          rbind,
          Map(function(x) {x$params}, private$children),
          c()
        )
        childParamsMean <- colMeans(childParams)
      }

      # Construct prior for node mean
      priorParamsMean <- (parentParams + numOfChildren*childParamsMean)/(numOfChildren + 1)
      priorParamsCov <- drift / (numOfChildren + 1)

      # Posterior for node mean
      if (numOfData == 0) {
        postParamsMean <- priorParamsMean
        postParamsCov <- priorParamsCov
      } else {
        tmpCov <- priorParamsCov + (1/numOfData) * self$sigma
        invTmpCov <- chol2inv(chol(tmpCov))
        postParamsCov <- priorParamsCov - priorParamsCov%*%invTmpCov%*%priorParamsCov
        postParamsMean <- (priorParamsMean%*%chol2inv(chol(priorParamsCov)) +
                             numOfData*dataMean%*%chol2inv(chol(self$sigma))) %*% postParamsCov
      }

      self$params <- rmvnorm(n = 1,
                             mean = postParamsMean,
                             sigma = postParamsCov)

      # Posterior for node covarirance
      if (numOfData > 0) {
        postSigmaDof <- priorSigmaDof + numOfData
        tmpVec <- nodeData - repmat(self$params, numOfData, 1)
        postSigmaScale <- priorSigmaScale + t(tmpVec) %*% tmpVec
        self$sigma <- riwish(v = postSigmaDof, S = postSigmaScale)
      }
      invisible(self)
    },

    ResampleHyperParams = function() {

      if(!is.null(private$parent)) {
        stop("Can only update hypers from root!")
      }
      # sample drift
      ComputeDriftLlh <- function(drift) {

        drift <- diag(drift, length(drift))

        if (sum(drift < 0) > 0) {
          return (-Inf)
        }
        Descend <- function(root) {
          llh <- 0
          children <- root$GetChildren()
          for (i in seq_along(children)) {
            child <- children[[i]]
            llh <- llh + dmvnorm(child$params, root$params, drift, log = TRUE)
            llh <- llh + Descend(child)
          }
          return(llh)
        }
        return(
          Descend(self) +
            dmvnorm(self$params, private$initMean, drift, log = TRUE) +
            diwish(drift, private$priorDriftDof, private$priorDriftScale))
      }
      tmp <- diag(private$drift)
      tmp <- SliceSampler(tmp, ComputeDriftLlh, stepOut = F)
      private$drift <- diag(tmp, length(tmp))
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