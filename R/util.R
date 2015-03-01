#' Convert beta random variables to edges in stick-breaking process
#' @param sticks A beta random variable vector
#' @return A vector of edges
#' @export
SticksToEdges <- function(sticks){
  return(1.0 - cumprod(1.0 - sticks))
}


#' Trim leading or trailing zeros in a vector
#' @param x A numerical vector
#' @param trim A string with "f" representing trim from front and "b" to
#' trim from back. Default is "fb", trim zeros from both front and back of the
#' array.
#' @return The result of trimming the input
#' @export
TrimZeros <- function(x, trim = "fb") {
  nonZeroIds <- which(x != 0)
  if (length(nonZeroIds) == 0) {
    return(c())
  }

  if (trim == "fb") {
    return(x[min(nonZeroIds):max(nonZeroIds)])
  } else if (trim == "f") {
    return(x[min(nonZeroIds):length(x)])
  } else {
    return(x[1:max(nonZeroIds)])
  }
}


ConvertFunctionNameToVariableName <- function(f) {
  res <- unlist(strsplit(gsub("(.)([[:upper:]])", "\\1 \\2", f), " "))
  res[2] <- paste(tolower(substring(res[2], 1, 1)), substring(res[2], 2), sep = "" )
  return(paste(res[-1], collapse = "" ))
}

# Doesn't work
GetHyperParamsTemplete <- function(f) {
  v <- ConvertFunctionNameToVariableName(f)
  function() {
    print(v)
    if (is.null(self$GetParent())) {
      return(eval(paste("private$", v, sep ="")))
    } else {
      return(do.call(paste("self$GetParent()$", f, sep = ""  ), list()))
    }
  }
}