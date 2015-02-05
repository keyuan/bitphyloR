#' R6 class for TSSB. TSSB is the basic object of tree structured stick-breaking
#' process
#'
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @keywords data


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
      self$root <- list(
        node     = rootNode,
        main     = if (self$minDepth == 0) rbeta(1, 1, self$dpAlpha) else 0,
        sticks   = c(),
        children = list())
      rootNode$tssb <- self

      # Draw data assignments and a tree
      for (n in 1:nrow(self$data)) {
        res <- self$FindNode(runif(1))
        self$assignments <- c(self$assignments, res$node)
        res$node$AddDatum(n)
        self$root <- res$root
      }
    },

    FindNode = function(u) {
      descend <- function(root, u, path = c(), depth = 0) {
        if (depth >= self$maxDepth) {
          warning("Reached maximum depth")
          return(list(node = root$node, path = path, root = root))
        } else if (u < root$main) {
          return(list(node = root$node, path = path, root = root))
        } else {
          # Rescale the uniform variate to the remaining interval.
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
                                  sticks = c(),
                                  children = list()
                                  )
                                )
                              )
          }
          edges <- SticksToEdges(root$sticks)
          index <- sum(u > edges) + 1
          edges <- c(0, edges)
          tempu <- u
          u     <- (u - edges[index]) / (edges[index+1] - edges[index])
          res <- descend(root$children[[index]], u, depth+1)
          node <- res$node
          path <- res$path
          path <- c(path,index)
          root$children[[index]] <- res$root
          return(list(node = node, path = path, root = root))
          }
      }
      return(descend(self$root, u))
    },

    GetMixture = function() {
      descend <- function(root, mass) {
        weight <- mass * root$main
        node <- root$node
        edges <- SticksToEdges(root$sticks)
        weights <- diff(c(0, edges))

        if (length(root$children) < 1) {
          return(list(node = node, weight = weight))
        } else {
          for (i in 1:length(root$children)) {
            child <- root$children[[i]]
            res <- descend(child, mass*(1.0-root$main)*weights[i])
            node <- c(node, res$node)
            weight <- c(weight, res$weight)
          }
          return(list(node = node, weight = weight))
        }
      }
      return(descend(self$root, 1.0))
    },

    ConvertTssbToIgraph = function() {
      edges <- SticksToEdges(self$root$sticks)
      weights <- diff(c(0, edges))
      g <- graph.empty(directed = TRUE)
      g <- g + vertex(name = "X",
                      size = self$root$node$GetNumOfLocalData())

      descend <- function(root, name, mass, g) {
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
            tmp <- descend(child,
                           childName,
                           mass*weights[i] * (1.0 - child$main),
                           g)
            g <- tmp$g
            total = total + childMass + tmp$total
          }
          return(list(total=total, g=g))
        }
      }
    res = descend(self$root, "X", 1, g)
    return(res)
    },

    CullTree = function() {}

  )
)


