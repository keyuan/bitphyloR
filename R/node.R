library(R6)

Node <- R6Class(
  classname = "Node",

  public = list(
    # fields
    dataIds = list(),
    tssb = emptyenv(),

    # methods
    initialize = function(parent = emptyenv(), tssb = emptyenv()) {
      if (is.environment(parent) &&
            !identical(parent, emptyenv()) &&
            class(parent)[1] == "Node") {
        private$parent <- parent
        parent$AddChild(self)
      } else {
        private$parent <- parent
      }
      if (is.environment(tssb) ) self$tssb <- tssb
      },


    GetChildren = function() {
      private$children
    },

    GetParent = function() {
      private$parent
    },

    AddChild = function(child) {
      if (!missing(child)) {
        private$children <- c(private$children, child)
      }
      invisible(self)
    },

    RemoveChild = function(child) {
      if (!missing(child)) {
        private$children <- Filter(
          Negate(
            function(x) identical(child, x)
            ),
          private$children)
      }
      invisible(self)
    },

    Kill = function() {
      if (!identical(private$parent, emptyenv())) {
        private$parent$RemoveChild(self)
      }

      private$children = NULL
      private$parent = NULL
    },

    Spawn = function() {
      Node$new(parent = self, tssb = self$tssb)
    }

    HasData = function() {
      if (length(self$dataIds)>0) {
        return(TRUE)
      } else {
        Reduce()
      }
      }
    }
#     def has_data(self):
#       if len(self.data):
#       return True
#     else:
#       for child in self._children:
#       if child.has_data():
#       return True
#     return False

    ), # end of public

  private = list(
    # fields
    children = list(),
    parent = emptyenv()
    )
  )

n0 <- Node$new()
n1 <- Node$new(parent = n0)
n2 <- Node$new(parent = n0)
