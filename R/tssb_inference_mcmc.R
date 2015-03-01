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

      for (n in 1:self$numOfData) {
        ancestors <- self$assignments[n]$GetAncestors()
        current <-  self$root
        indices <- c()
        for (anc in seq_along(ancestors)) {
          index = which(
            unlist(
              Map(function(x) identical(x,anc), current$children)
            )
          )
          current = current$children[index]
          indices = c(indices, index)
        }
        maxU <- 1
        minU <- 0
        llhS = log(runif(1)) + self$assignments[n]$LogProb(self$data[n,])
        while (T) {
          newU <- (maxU - minU)*runif(1) + minU
          res <- self$FindNode(newU)
          newNode <- res$node
          newPath <- res$path
          newLlh <- newNode$LogProb(self$data[n,])
          if (newLlh > llhS) {
            if (!identical(newNode != self$assignments[n])) {
              self$assignments[n]$RemoveDatum(n)
              newNode$AddDatum(n)
              self.assignments[n] = newNode
              break
            }
          } else if (abs(maxU - minU) < .Machine$double.eps) {
            warning("Slice sampler shrank down. Keep current state.")
            break
          } else {
            res <- ComparePath(indices, newPath)
            if (res < 0) {
              minU <- newU
            } else if (res >0) {
              maxU <- newU
            } else {
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
          root$sticks[i] <- rbeta(1, postAlpha, postBeta)
          dataDown <- dataDown + childData
        }

        dataHere <- root$node$GetNumOfLocalData()
        postAlpha <- 1 + dataHere
        postBeta <- (self$dpLambda^depth) * self$dpAlpha + dataDown
        root$main = if (self$minDepth <= depth) rbeta(1, postAlpha, postBeta) else 0.0
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
              root$sticks <- c(root$sticks, rbeta(1, 1, self$dpGamma))
              root$children <- c(root$children,
                                 list(
                                   list(
                                     node = root$node$Spawn(),
                                     main = if (self$minDepth <= depth +1 ) {
                                       rbeta(1, 1, (self$dpLambda^(depth+1))*self$dpAlpha)
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
          child = root$children[k]
          newChildren <- c(newChildren, child)
          root$children[[k]] <- Descend(root$children[[k]], depth + 1)
        }

        lapply(Filter(function(x) ! x %in% newOrder, seq_along(root$sticks)),
               function(k) {root$children[[k]]$node$Kill()
                            root$children[[k]] <- list()})

        root$children = newChildren
        root$sticks   = rep(0, length(newChildren))

        return(root)
      }
      self$root <- Descend(self$root)
      self$ResampleSticks()
      invisible(self)
    }

  )
)





