#'@include tssb.R
NULL

#' R6 class for inference via MCMC for TSSB.
#'
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @keywords data

TssbMCMC <- R6::R6Class(
  classname = "TssbMCMC",
  inherit = TSSB,
  public = list(
    initialize = function(...) {
      super$initialize(...)
    },

    ResampleAssignments = function() {
      lengths <- c()

      for (n in seq_len(self$numOfData)) {
        ancestors <- self$assignments[[n]]$GetAncestors()
        current <-  self$root
        indices <- c()
        for (anc in seq_along(ancestors[-1])) {
          index = which(
            unlist(
              Map(function(x) identical(x$node, ancestors[[anc+1]]), current$children)
            )
          )
          current = current$children[[index]]
          indices = c(indices, index)
        }
        maxU <- 1
        minU <- 0
        llhS = log(runif(1)) + self$assignments[[n]]$GetLogProb(self$data[n,])
        while (T) {
          newU <- (maxU - minU)*runif(1) + minU
          res <- self$FindNode(newU)
          newNode <- res$node
          newPath <- res$path
          self$root <- res$root
          newLlh <- newNode$GetLogProb(self$data[n,])
          if (newLlh > llhS) {
            if (!identical(newNode, self$assignments[[n]])) {
              self$assignments[[n]]$RemoveDatum(n)
              newNode$AddDatum(n)
              self$assignments[[n]] <- newNode
            }
            break
          } else if (abs(maxU - minU) < .Machine$double.eps) {
            warning("Slice sampler shrank down. Keep current state.")
            break
          } else {
            pathCmp <- ComparePath(indices, newPath)
            if (pathCmp < 0) {
              minU <- newU
            } else if (pathCmp >0) {
              maxU <- newU
            } else {
              cat("id: ", n, " indices: ",indices,  " newPath: ", newPath,
                  " llhS: ", llhS, " llhNew: ", newLlh,
                  " newU: ", newU, " minU: ", minU, " maxU: ", maxU)
              stop("Slice sampler weirdness")
            }
          }
        }
      }
      invisible(self)
    },

    ResampleSticks = function() {
      Descend <- function(root, depth=0) {
        dataDown <- 0
        indices <- seq_along(root$children)

        for (i in sort(indices, decreasing = T)) {
          res <- Descend(root$children[[i]], depth+1)
          childData <- res$nodeData
          root$children[[i]] <- res$root
          postAlpha <- 1 + childData
          postBeta <- self$dpGamma + dataDown
          root$sticks[i] <- BoundBeta(1, postAlpha, postBeta)
          dataDown <- dataDown + childData
        }

        dataHere <- root$node$GetNumOfLocalData()
        postAlpha <- 1 + dataHere
        postBeta <- (self$dpLambda^depth) * self$dpAlpha + dataDown
        root$main = if (self$minDepth <= depth) BoundBeta(1, postAlpha, postBeta) else 0.0
        return(list(nodeData = dataHere + dataDown, root = root))
      }

      self$root = Descend(self$root)$root
      invisible(self)
    },

    ResampleStickOrders = function() {
      Descend <- function(root, depth = 0) {
        if (length(root$children)==0) {
          return(root)
        }
        newOrder <- c()
        represented <- Filter(function(i) root$children[[i]]$node$HasData(),
                              seq_along(root$children))
        allWeights <- diff(c(0,SticksToEdges(root$sticks)))

        while (length(represented)!=0) {
          u <- runif(1)
          while (TRUE) {
            subIndices <- Filter(function(x) ! x %in% newOrder, seq_along(root$sticks))
            subWeights <- c(allWeights[subIndices], 1-sum(allWeights))
            subWeights <- subWeights/sum(subWeights)
            index <- sum(u > cumsum(subWeights))
            if (index == length(subIndices)) {
              root$sticks <- c(root$sticks, BoundBeta(1, 1, self$dpGamma))
              root$children <- c(root$children,
                                 list(
                                   list(
                                     node = root$node$Spawn(),
                                     main = if (self$minDepth <= depth +1 ) {
                                       BoundBeta(1, 1, (self$dpLambda^(depth+1))*self$dpAlpha)
                                     } else {0},
                                     sticks = c(),
                                     children = list())))
              allWeights <- diff(c(0, SticksToEdges(root$sticks)))
            } else {
              index <- subIndices[index+1]
              break
            }
          }
          newOrder <- c(newOrder, index)
          represented <- represented[represented != index]
        }

        newChildren <- c()
        for (k in newOrder) {
          child <- root$children[k]
          newChildren <- c(newChildren, child)
          root$children[[k]] <- Descend(root$children[[k]], depth + 1)
        }

        lapply(Filter(function(x) ! x %in% newOrder, seq_along(root$sticks)),
               function(k) {root$children[[k]]$node$Kill()
                            root$children[[k]] <- list()})

        root$children <- newChildren
        root$sticks <- rep(0, length(newChildren))

        return(root)
      }
      self$root <- Descend(self$root)
      self$ResampleSticks()
      invisible(self)
    },

    ResampleNodeParameters = function() {
      Descend <- function(root) {
        Map(Descend, root$children)
        root$node$ResampleParams()
        return(root=root)
      }
      self$root <- Descend(self$root)
      invisible(self)
    },

    ResampleHypers = function(sampleDpAlpha = T,
                              sampleDpLambda = T,
                              sampleDpGamma = T) {

      ComputeDpMainLlh <- function(params,
                                   sample = c("alpha", "lambda"),
                                   fixParams = NULL) {
        if (sample == "alpha") {
          dpAlpha <- params
          dpLambda <- fixParams
        } else if (sample == "lambda") {
          dpAlpha <- fixParams
          dpLambda <- params
        } else {
          stop("sample has to be either 'alpha' or 'lambda' ")
        }

        if (dpAlpha < self$minDpAlpha || dpAlpha > self$maxDpAlpha) {
          return(-Inf)
        }
        if (dpLambda < self$minDpLambda || dpLambda > self$maxDpLambda) {
          return(-Inf)
        }
        Descend <- function(root, depth = 0) {
          llh <- if (self$minDepth <= depth) {
            dbeta(root$main, 1, dpLambda^depth*dpAlpha, log = TRUE)
          } else {0}

          for (i in seq_along(root$children)) {
            res <- Descend(root$children[[i]], depth + 1)
            llh <- llh + res$llh
            root$children[[i]] = res$root
          }
          return(list(llh = llh, root = root))
        }
        return(Descend(self$root)$llh)
      }

      ComputeDpBranchLlh <- function(dpGamma) {
        if (dpGamma < self$minDpGamma || dpGamma > self$maxDpGamma) {
          return(-Inf)
        }
        Descend <- function(root) {
          llh <- 0
          for (i in seq_along(root$children)) {
            llh <- llh + dbeta(root$sticks[i], 1, dpGamma, log = TRUE)
            res <- Descend(root$children[[i]])
            llh <- llh + res$llh
            root$children[[i]] <- res$root
          }
          return(list(llh = llh, root = root))
        }
        return(Descend(self$root)$llh)
      }

      if (sampleDpAlpha) {
        self$dpAlpha <- SliceSampler(self$dpAlpha,
                                     ComputeDpMainLlh,
                                     sample = "alpha",
                                     fixParams = self$dpLambda)
      }
      if (sampleDpLambda) {
        self$dpLambda <- SliceSampler(self$dpLambda,
                                      ComputeDpMainLlh,
                                      sample = "lambda",
                                      fixParams = self$dpAlpha)

      }
      if (sampleDpGamma) {
        self$dpGamma <- SliceSampler(self$dpGamma, ComputeDpBranchLlh)
      }
      invisible(self)
    }
  )
)





