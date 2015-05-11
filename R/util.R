
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

#' Sample from a user supplied distribution with slice sampler
#' @param initX initial parameter
#' @param logprob a user supplied log-distribution function
#' @param sigma step size for slice construction
#' @param setpOut use setp out procedure or not
#' @param maxStepsOut maximum step out number
#' @param compwise use component wise update or not
#' @param ... additional args for logprob
#' @return a sample of parameter
#' @export
SliceSampler <- function(initX,
                         logprob,
                         sigma = 1.0,
                         stepOut = T,
                         maxStepsOut = 1000,
                         compwise = F, ...) {

  DirectionSlice <- function(direction, initX) {

    DirLogProb <- function(z) {
      return(logprob(direction*z + initX, ...))
    }

    llhS <- log(runif(1)) + DirLogProb(0)
    upper <- sigma * runif(1)
    lower <- upper - sigma
    lowerStepsOut <- 0
    upperStepsOut <- 0
    if (stepOut) {

      while (DirLogProb(lower) > llhS && lowerStepsOut < maxStepsOut) {
        lowerStepsOut <- lowerStepsOut + 1
        lower <- lower - sigma
      }
      while (DirLogProb(upper) > llhS && upperStepsOut < maxStepsOut) {
        upperStepsOut <- upperStepsOut + 1
        upper <- upper + sigma
      }
    }

    stepIn <- 0
    while(T) {
      stepIn <- stepIn + 1
      newZ <- (upper - lower) * runif(1) + lower
      llhNew <- DirLogProb(newZ)
      if (is.nan(llhNew)) {
        cat("%f, %f, %f", newZ, initX, direction*newZ + initX, logprob(initX, ...))
        stop("Slice sampler got a NaN")
      }
      if (llhNew > llhS) {
        break
      } else if (newZ < 0) {
        lower <- newZ
      } else if (newZ > 0) {
        upper <- newZ
      } else {
        stop("Slice sampler shrank to zero")
      }
    }
    return(direction*newZ + initX)
  }

  dims <- length(initX)
  if (compwise) {
    ordering <- sample(dims)
    newX <- Reduce(
      rbind,
      Map(
        function(d) {
          direction = array(0, dim = dims)
          direction[d] = 1
          return(direction*DirectionSlice(direction, initX))
        },
        ordering),
      c())
    return(colSums(newX))
  } else {
    direction <- rnorm(dims)
    direction <- direction / sqrt(sum(direction^2))
    return(DirectionSlice(direction, initX))
  }
}


repmat <- function(x, n, m){
  x <- as.matrix(x)
  return(kronecker(matrix(1,n,m),x))
}

ConvertFunctionNameToVariableName <- function(f) {
  res <- unlist(strsplit(gsub("(.)([[:upper:]])", "\\1 \\2", f), " "))
  res[2] <- paste(tolower(substring(res[2], 1, 1)), substring(res[2], 2), sep = "" )
  return(paste(res[-1], collapse = "" ))
}

cmp <- function(x, y) {
  if (x > y) {
    return(1)
  } else if (x < y) {
    return(-1)
  } else {
    return(0)
  }
}

ComparePath <- function(x, y) {
  if (is.null(x) && is.null(y)) {
    return(0)
  } else if (is.null(x)) {
    return(1)
  } else if (is.null(y)) {
    return(-1)
  }
  s1 <- paste(as.character(x), collapse = "")
  s2 <- paste(as.character(y), collapse = "")
  return(cmp(s2, s1))
}


