
# TSSB  -------------------------------------------------------------------

#' R6 class for TSSB. TSSB is the basic object of tree structured stick-breaking
#' process
#'
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @keywords data
#' @field data
#' @field dpAlpha
#' @field dpGamma
#' @field dpLambda
#' @field minDepth
#' @field maxDepth
TSSB <- R6Class(
  classname = "TSSB",
  public = list(
    # Felids ------------------------------------------------------------------
    data = NULL,
    dpAlpha = NA,
    dpGamma = NA,
    dpLambda = NA,
    minDepth = NA,
    maxDepth = NA,
    root = list(),
    assignments = list(),
    numOfData = NA,

    # Methods -----------------------------------------------------------------
    initialize = function(rootNode = emptyenv(),
                          dpAlpha = 1,
                          dpGamma = 1,
                          dpLambda = 1,
                          data = NULL,
                          minDepth = 0,
                          maxDepth = 25
                          ) {
      if (!is.environment(rootNode) ||
            identical(rootNode, emptyenv()) ||
            class(rootNode)[1] != "Node") {
        stop("Root node must be specified")
      }

      self$data    <- data
      self$dpAlpha <- dpAlpha
      self$dpGamma <- dpGamma
      self$minDepth <- minDepth
      self$maxDepth <- maxDepth
      self$dpLambda <- dpLambda
      self$numOfData <- if (is.null(data)) 0 else nrow(self$data)
      self$root <- list(
        node     = rootNode,
        main     = if (self$minDepth == 0) rbeta(1, 1, self$dpAlpha) else 0,
        sticks   = c(),
        children = list())
      rootNode$tssb <- self

      # Draw data assignments and a tree
      uSamples <- runif(self$numOfData)
      for (n in 1:self$numOfData) {
        res <- self$FindNode(uSamples[n])
        self$assignments <- c(self$assignments, res$node)
        res$node$AddDatum(n)
        self$root <- res$root
      }
    },

    FindNode = function(u) {
      #' Top-down transverse the tree
      Descend <- function(root, u, path = c(), depth = 0) {
        if (depth >= self$maxDepth) {
          warning("Reached maximum depth")
          return(list(node = root$node, path = path, root = root))
        } else if (u < root$main) {
          return(list(node = root$node, path = path, root = root))
        } else {
          #' Rescale the uniform variate to the remaining interval.
          u <- (u - root$main) / (1 - root$main)
          while (length(root$children) == 0
                 || u > (1 - prod(1 - root$sticks))
                 ) {
            root$sticks <- c(root$sticks, rbeta(1, 1, self$dpGamma))
            root$children <- c(root$children,
                              list(
                                list(
                                  node = root$node$Spawn(),
                                  main = if (self$minDepth <= (depth+1)) {
                                    rbeta(1, 1, (self$dpLambda^(depth+1))*self$dpAlpha)
                                    } else {0},
                                  sticks = NULL,
                                  children = list()
                                  )
                                )
                              )
          }
          edges <- c(0, SticksToEdges(root$sticks))
          index <- sum(u > edges)
          u     <- (u - edges[index]) / (edges[index+1] - edges[index])
          res <- Descend(root$children[[index]], u, path, depth+1)
          node <- res$node
          path <- c(index, res$path)
          root$children[[index]] <- res$root
          return(list(node = node, path = path, root = root))
          }
      }
      return(Descend(self$root, u))
    },


    GetMixture = function() {
      Descend <- function(root, mass) {
        weight <- mass * root$main
        node <- root$node
        edges <- SticksToEdges(root$sticks)
        weights <- diff(c(0, edges))

        for (i in seq_along(root$children)) {
          child <- root$children[[i]]
          res <- Descend(child, mass*(1.0-root$main)*weights[i])
          node <- c(node, res$node)
          weight <- c(weight, res$weight)
        }
        return(list(node = node, weight = weight))
      }
      return(Descend(self$root, 1.0))
    },

    ConvertTssbToIgraph = function() {
      edges <- SticksToEdges(self$root$sticks)
      weights <- diff(c(0, edges))
      g <- graph.empty(directed = TRUE)
      g <- g + vertex(name = "X",
                      size = self$root$node$GetNumOfLocalData())

      Descend <- function(root, name, mass, g) {
        if (length(root$sticks) < 1){
          return(list(total = mass, g=g))
        } else {
          total <- 0
          edges <- SticksToEdges(root$sticks)
          weights <- diff(c(0, edges))

          for (i in 1:length(root$children)) {
            child <- root$children[[i]]
            childName <- paste(name, i, sep = "-")
            childMass <- mass * weights[i] * child$main
            g <- g + vertex(name = childName,
                            size = child$node$GetNumOfLocalData())
            g <- g + edge(name, childName,
                          Value = child$node$GetNumOfLocalData())
            tmp <- Descend(child,
                           childName,
                           mass*weights[i]*(1.0 - child$main),
                           g)
            g <- tmp$g
            total = total + childMass + tmp$total
          }
          return(list(total=total, g=g))
        }
      }
    res = Descend(self$root, "X", 1, g)
    return(res)
    },

    CullTree = function() {
      Descend <- function(root) {
        res <- unlist(Map(Descend, root$children), recursive = F, use.names = F)

        if (length(res) == 0) {
          return(list(counts = root$node$GetNumOfLocalData(), root = root))
        }
        counts <- unlist(res[seq(1, length(res), by = 2)])
        root$children <- res[seq(2, length(res), by = 2)]
        keep <- which(counts != 0)

        sapply(root$children[which(counts == 0)], function(x) x$node$Kill())

        if (length(keep) == 0) {
          root$sticks <- NULL
          root$children <- list()
        } else {
          root$sticks <- root$sticks[keep]
          root$children <- root$children[keep]
        }
        return(list(counts = sum(counts) + root$node$GetNumOfLocalData(),
                    root = root))
      }
      res <- Descend(self$root)
      self$root <- res$root
      invisible(self)
    }
  )
)

# MCMC --------------------------------------------------------------------

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
  inherit = TSSB,
  public = list(
    ResampleAssignments = function() {
      epsilon <- 2.2204460492503131e-16
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
          } else if (abs(maxU - minU) < epsilon) {
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
      invisible(self)
    }

    )
)


# Batch VB ----------------------------------------------------------------


#' R6 class for inference via Batch variational Bayes for TSSB.
#'
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @keywords data

TssbBatchVB <- R6Class(
  classname = "TssbBatchVB",
  inherit = TSSB
)


# Stochastic VB -----------------------------------------------------------


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

