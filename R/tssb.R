TSSB <- R6Class(
  classname = "TSSB",
  public = list(
    # Felids ------------------------------------------------------------------
    data = NULL,
    dpAlpha = NA,
    dpGamma = NA,
    minDepth = NA,
    maxDepth = NA,
    root = list(),

    # Methods -----------------------------------------------------------------
    initialize = function(rootNode = emptyenv(),
                          dpAlpha = 1,
                          dpGamma = 1,
                          data = NULL,
                          minDepth = 0,
                          maxDepth = 15
                          ) {
      if (!is.environment(rootNode) ||
            identical(rootNode, emptyenv()) ||
            class(rootNode)[1] != "Node") {
        stop("Root node must be specified")
      } else {
        self$data    = data
        self$dpAlpha = dpAlpha
        self$dpGamma = dpGamma
        self$minDepth = minDepth
        self$root = list(
          "node"     = rootNode,
          "main"     = if (self$minDepth == 0) rbeta(1, 1.0, self$dpAlpha) else 0,
          "sticks"   = c(),
          "chdren" = list())
        rootNode$tssb = self
      }
    },

    FindNode = function(u) {
      desend = function(root, u, depth =0 ) {
        if (depth >= self$maxDepth) {
          warning("Reached maximum depth")
          return(list(node = root$node, path = c()))
        } else if (u < root$main) {
          return(list(node = root$node, path = c()))
        } else {
          # Rescale the uniform variate to the remaining interval.
          u = (u - root$main) / (1 - root$main)
          while (length(root$children == 0)
                 || u > (1 - prod(1 - root$sticks))
                 ) {
            root$sticks = c(root$sticks, rbeta(1, self$dpGamma))
          }

        }

      }
    }



    def find_node(self, u):
      def descend(root, u, depth=0):
      if depth >= self.max_depth:
      #print >>sys.stderr, "WARNING: Reached maximum depth."
      return (root['node'], [])
    elif u < root['main']:
      return (root['node'], [])
    else:
      # Rescale the uniform variate to the remaining interval.
      u = (u - root['main']) / (1.0 - root['main'])

    # Perhaps break sticks out appropriately.
    while not root['children'] or (1.0 - prod(1.0 - root['sticks'])) < u:
      root['sticks'] = vstack([ root['sticks'], boundbeta(1, self.dp_gamma) ])
    root['children'].append({ 'node'     : root['node'].spawn(),
                              'main'     : boundbeta(1.0, (self.alpha_decay**(depth+1))*self.dp_alpha) if self.min_depth <= (depth+1) else 0.0,
                              'sticks'   : empty((0,1)),
                              'children' : [] })

    edges = 1.0 - cumprod(1.0 - root['sticks'])
    index = sum(u > edges)
    edges = hstack([0.0, edges])
    u     = (u - edges[index]) / (edges[index+1] - edges[index])

    (node, path) = descend(root['children'][index], u, depth+1)

    path.insert(0, index)

    return (node, path)
    return descend(self.root, u)

#     root_node.tssb = self
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
