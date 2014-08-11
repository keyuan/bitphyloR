Node <- R6Class(
  classname <- "Node",

  # Public ====================================================================#
  public <- list(
    ## Fields ----------------------------------------------------------------##
    data <- list(),
    tssb <- NA,
    children <- list(),
    parent <- NA,
    ## -----------------------------------------------------------------------##

    ## Methods ---------------------------------------------------------------##
    initialize = function(parent = NA, tssb = NA) {
      self$parent <- parent
      self$tssb <- tssb
      },

    AddChild <- function(child) {
      if (!missing(child)) self$children <- c(self$children, child)
      invisible(self)
    }
    ## end of methods --------------------------------------------------------##

    )
  # End of public =============================================================#
  )


