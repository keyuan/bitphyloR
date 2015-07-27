#'TSSB is a R6 object of tree-structured stick breaking
#' process
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
#' @field root
#' @field assignments
#' @field numOfData
#' @field minDpAlpha
#' @field maxDpAlpha
#' @field minDpLambda
#' @field minDpLambda
#' @field minDpGamma
#' @field maxDpGamma
#' @method new TSSB constuctor
#' @method FindNode Find node in tree, generate new node as needed
#' @method GetMixture Compute the mixing weights and get the corresponding nodes
#' @method ConvertTssbToIgraph Convert the tssb root list to an igraph object
#' @method CullTree Remove empty leaf nodes
#' @method GetLogMarginalDataLikelihood Compute the log marginal data likelihood
TSSB <- R6::R6Class(
  classname = "TSSB",
  public = list(
    data = NULL,
    # Felids ------------------------------------------------------------------
    dpAlpha = NA,
    dpGamma = NA,
    dpLambda = NA,
    minDepth = NA,
    maxDepth = NA,
    root = NULL,
    assignments = NULL,
    numOfData = NA,
    minDpAlpha = NA,
    maxDpAlpha = NA,
    minDpLambda = NA,
    maxDpLambda = NA,
    minDpGamma = NA,
    maxDpGamma = NA,

    # Methods -----------------------------------------------------------------
    initialize = function(rootNode = emptyenv(),
                          dpAlpha = 1,
                          dpGamma = 1,
                          dpLambda = 1,
                          data = NULL,
                          minDepth = 0,
                          maxDepth = 25,
                          minDpAlpha = 1e-6, maxDpAlpha = 50,
                          minDpLambda = 1e-6, maxDpLambda = 1,
                          minDpGamma = 1e-6, maxDpGamma = 50
                          ) {
      if (!is.environment(rootNode) ||
            identical(rootNode, emptyenv()) ||
            !"Node" %in% class(rootNode)) {
        stop("Root node must be specified")
      }

      self$data    <- data
      self$dpAlpha <- dpAlpha
      self$dpGamma <- dpGamma
      self$minDepth <- minDepth
      self$maxDepth <- maxDepth
      self$dpLambda <- dpLambda
      self$minDpAlpha <- minDpAlpha
      self$maxDpAlpha <- maxDepth
      self$minDpLambda <- minDpLambda
      self$maxDpLambda <- maxDpLambda
      self$minDpGamma <- minDpGamma
      self$maxDpGamma <- maxDpGamma
      self$numOfData <- if (is.null(data)) 0 else nrow(self$data)
      self$root <- list(
        node     = rootNode,
        main     = if (self$minDepth == 0) BoundBeta(1, 1, self$dpAlpha) else 0,
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
            root$sticks <- c(root$sticks, BoundBeta(1, 1, self$dpGamma))
            root$children <- c(root$children,
                              list(
                                list(
                                  node = root$node$Spawn(),
                                  main = if (self$minDepth <= (depth+1)) {
                                    BoundBeta(1, 1, (self$dpLambda^(depth+1))*self$dpAlpha)
                                    } else {0},
                                  sticks = NULL,
                                  children = NULL
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
      g <- igraph::graph.empty(directed = TRUE)
      g <- g + igraph::vertex(name = "X",
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
            g <- g + igraph::vertex(name = childName,
                                    size = child$node$GetNumOfLocalData())
            g <- g + igraph::edge(name, childName,
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
          root$children <- NULL
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
    },

    GetLogMarginalDataLikelihood = function() {
      res <- self$GetMixture()
      weights <- res$weight
      nodes <- res$node

      ll <- Reduce(
        sum,
        Map(function(i) {
          if (length(weights) == 1) {
            node <- nodes
          } else {
            node <- nodes[[i]]
          }
          if (node$GetNumOfLocalData()) {
            node$GetNumOfLocalData()*weights[i] + node$GetNodeLogProb()
          } else {
            0
          }
        },
        seq_along(weights)
        ),
        0
      )
      return(list(ll=ll, ww = weights, nn = nodes))
    }
  )
)



