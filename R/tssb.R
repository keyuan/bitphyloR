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
    }

    
#     
#     def __init__(self, dp_alpha=1.0, dp_gamma=1.0, root_node=None, data=None,
#                  min_depth=0, max_depth=15, alpha_decay=1.0):
#       if root_node is None:
#       raise Exception("Root node must be specified.")
#     
#     self.min_depth   = min_depth
#     self.max_depth   = max_depth
#     self.dp_alpha    = dp_alpha
#     self.dp_gamma    = dp_gamma
#     self.alpha_decay = alpha_decay
#     self.data        = data
#     self.num_data    = 0 if data is None else data.shape[0]
#     self.root        = { 'node'     : root_node,
#                          'main'     : boundbeta(1.0, dp_alpha) if self.min_depth == 0 else 0.0,
#                          'sticks'   : empty((0,1)),
#                          'children' : [] }
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
