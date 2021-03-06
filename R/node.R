#' Node is a R6 object of each cluster in TSSB
#'
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @keywords Node class
#' @field dataIds: data IDs
#' @field tssb: A TSSB object
#' @field children: descentant node objects
#' @field parent: parent node object
#' @method new
#' @method GetChildren
#' @method GetParent
#' @method SetParent
#' @method AddChild
#' @method RemoveChild
#' @method Kill
#' @method Spawn
#' @method HasData
#' @method AddDatum
#' @method RemoveDatum
#' @method GetNumOfLocalData
#' @method GetNumOfSubTreeData
#' @method GetData
#' @method GetAncestors

Node <- R6::R6Class(
  classname = "Node",

  public = list(
    # Fields ------------------------------------------------------------------
    dataIds = c(),
    tssb = NULL,

    # Methods -----------------------------------------------------------------
    initialize = function(parent = NULL, tssb = NULL) {
      if (!is.null(parent)) {
        private$parent <- parent
        parent$AddChild(self)
      } else {
        private$parent <- parent
      }
      if (!is.null(tssb) ) self$tssb <- tssb
      },


    GetChildren = function() {
      private$children
    },

    GetParent = function() {
      private$parent
    },

    SetParent = function(parent = NULL) {
      private$parent <- parent
    },

    AddChild = function(child) {
      if (!missing(child)) {
        child$SetParent(self)
        private$children <- c(private$children, child)
      }
      invisible(self)
    },

    RemoveChild = function(child) {
      if (!missing(child)) {
        private$children <- unlist(Filter(
          Negate(
            function(x) identical(child, x)
            ),
          private$children))
      }
      invisible(self)
    },

    Kill = function() {
      if (!is.null(private$parent))  {
        private$parent$RemoveChild(self)
      }
      private$children = NULL
      private$parent = NULL
    },

    Spawn = function() {
      return(get(class(self)[1])$new(parent = self, tssb = self$tssb))
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
      if (!missing(id) && !any(self$dataIds==id)) {
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
      if (dim(self$tssb$data)[2] > 1) {
        return(self$tssb$data[self$dataIds,])
      } else {
        return(as.matrix(self$tssb$data[self$dataIds,]))
      }

    },

    GetAncestors = function() {
      ancestors = c()
      if (is.null(private$parent)) {
        return(list(self))
      } else {
        ancestors <- c(private$parent$GetAncestors(), self)
        return(ancestors)
      }
    }

    ), # end of public

  private = list(
    # Fields
    children = NULL,
    parent = NULL
    )
  )




