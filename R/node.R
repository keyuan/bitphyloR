Node <- R6Class(
  classname = "Node",

  public = list(

    # fields
    dataIds = list(),
    tssb = NA,

    # methods
    initialize = function(parent = NA, tssb = NA) {
      if (!missing(parent)) private$parent <- parent
      if (!missing(tssb)) self$tssb <- tssb
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
    }
    ), # end of public

  private = list(
    # fields
    children = list(),
    parent = NA
    )
  )

n0 <- Node$new()
n1 <- Node$new(parent = n0)
n2 <- Node$new(parent = n0)
