source("R/node.R")

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
                          maxDepth = 15
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
        u <- runif(1)
        res <- self$FindNode(u)
        self$assignments <- c(self$assignments, res$node)
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
          root <- res$root
          return(list(node = node, path = path, root = root))
          }
      }
      return(descend(self$root, u))
    }

#
#     if False:
#       data_u           = rand(self.num_data)
#     self.assignments = []
#     for n in range(self.num_data):
#       (c, path) = self.find_node(data_u[n])
#     c.add_datum(n)
#     self.assignments.append(c)
#     else:
#       self.assignments = []
#     for n in range(self.num_data):
#       self.root['node'].add_datum(n)
#     self.assignments.append(self.root['node'])
  )
)

tssb <- TSSB$new(n0, data = matrix(rnorm(50),50,1))

