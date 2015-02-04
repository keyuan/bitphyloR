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
    root1 = list(),
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
      descend <- function(root, u, depth = 0) {
        if (depth >= self$maxDepth) {
          warning("Reached maximum depth")
          return(list(node = root$node, path = c(), root = root))
        } else if (u < root$main) {
          return(list(node = root$node, path = c(), root = root))
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
          edges <- 1 - cumprod(1 - root$sticks)
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
                      size = length(self$root$node$GetData()))

      descend <- function(root, name, mass, g) {
        total <- 0
        edges <- SticksToEdges(root$sticks)
        weights <- diff(c(0, edges))

        for (i in 1:length(root$children)) {
          child <- root$children[[i]]
          childName <- paste(name, i, sep = "-")
          childMass <- mass * weights[i] * child$main
          g <- g + vertex(name = childName,
                          size = length(child$node$GetData()))
          g <- g + edge(name, childName, Value = length(child$node$GetData()))
          tmp <- descend(child,
                         childName,
                         mass*weights[i] * (1.0 - child$main),
                         g)
          g <- tmp$g
          total = total + child_mass + tmp$total
        }
        return(list(total, g))
      }
      res = descend(self$root, "X", 1, g)
      return(res$g)
    }
  )
)

# tssb <- TSSB$new(n0, data = matrix(rnorm(50),10000,1))
# res <- tssb$GetMixture()

# def get_mixture(self):
#   def descend(root, mass):
#     weight  = [ mass * root['main'] ]
#     node    = [ root['node'] ]
#     edges   = sticks_to_edges(root['sticks'])
#     weights = diff(hstack([0.0, edges]))
#
#     for i, child in enumerate(root['children']):
#       (child_weights, child_nodes) = descend(child, mass*(1.0-root['main'])*weights[i])
#     weight.extend(child_weights)
#     node.extend(child_nodes)
#     return (weight, node)
#   return descend(self.root, 1.0)

def tssb2igraph(self):

  edges   = sticks_to_edges(self.root['sticks'])
  weights = diff(hstack([0.0, edges]))
  if len(weights) > 0:
    root_mass = weights[0] * self.root['main']
  else: #in case there is a problem with weights, as observed in some runs
    root_mass = -1

  g = Graph(directed = True)
  g.add_vertex(name = "X",
              params = " ".join(map(lambda x: "%.15f" %x,
                                     self.root['node'].params)),
              size = len(self.root['node'].get_data()),
              mass = root_mass,
              branch = self.root['node'].branch_length,
              members = " ".join(map(lambda x: "%s" %x,
                                      self.root['node'].data)) )
