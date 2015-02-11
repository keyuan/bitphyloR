#' R6 class for Node. Node is the basic object of each cluster in TSSB
#'
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @keywords data
#' @field dataIds: data IDs
#' @field tssb: A TSSB object
#' @field children: descentant node objects
#' @field parent: parent node object
#' @method initialize
Node <- R6Class(
  classname = "Node",

  public = list(
    # Fields ------------------------------------------------------------------
    dataIds = c(),
    tssb = emptyenv(),

    # Methods -----------------------------------------------------------------
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

    SetParent = function(parent = emptyenv()) {
      private$parent <- parent
    }

    AddChild = function(child) {
      if (!missing(child)) {
        child$SetParent(self)
        private$children <- c(private$children, child)
      }
      invisible(self)
    },

    RemoveChild = function(child) {
      if (!missing(child)) {
        child$SetParent()
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
        sapply(private$children, private$parent$AddChild)
      }
      private$children = NULL
      private$parent = NULL
    },

    Spawn = function() {
      return(Node$new(parent = self, tssb = self$tssb))
    },

    HasData = function() {
      if ( length(self$dataIds)>0 ) {
        return(TRUE)
      } else {
        return(
          sum(
            unlist(sapply(private$children, function(x) x$HasData()))
            )
          > 0)
      }
    },

    AddDatum = function(id) {
      if (!missing(id) && !id %in% self$dataIds) {
        self$dataIds <- c(self$dataIds, id)
      }
    },

    RemoveDatum = function(id) {
      if (!missing(id)) {
        if (!id %in% self$dataIds) {
          warning("id is not found in dataIds, nothing is removed")
        } else {
          self$dataIds <- self$dataIds[self$dataIds != id]
        }
      }
    },

    GetNumOfLocalData = function() {
      length(self$dataIds)
    },

    GetNumOfSubTreeData = function() {
        Reduce(
          sum,
          Map(function(x) x$GetNumOfSubTreeData(),
              private$children),
          length(self$dataIds)
        )
    },

    GetData = function() {
      self$tssb$data[self$dataIds,]
    },

    GetAncestors = function() {
      ancestors = c()
      if (identical(private$parent, emptyenv())) {
        return(self)
      } else {
        ancestors <- c(private$parent$GetAncestors(), self)
        return(ancestors)
      }

    }

    ), # end of public

  private = list(
    # Fields
    children = list(),
    parent = emptyenv()
    )
  )


