#' Fit a tssb with normal node with data via mcmc
#' @param data a data matrix or list
#' @param numOfMCMC mcmc iterations
#' @param burnIn burn-in number
#' @param dims dimension of data
#' @param tssb a use supplied TssbMCMC object
#' @param ... additional parameters to node class generator
#' @return A list of traces
#' @export

RunNormal <- function (data,
                       numOfMCMC = 1000,
                       burnIn = 0,
                       dims = ncol(data),
                       tssb = NULL,
                       ...) {

  if (is.list(data)) {
    numOfData <- length(data)
  } else if (is.matrix(data)) {
    numOfData <- nrow(data)
  } else {
    stop("data has to be list or matrix")
  }

  gg = vector(mode = "list", length = numOfMCMC)
  ll = array(0, dim = numOfMCMC)
  pp = array(0, dim = c(m, dims, numOfMCMC))
  ss = array(0, dim = c(m*dims, dims, numOfMCMC))
  dpa = array(0, dim = numOfMCMC)
  dpl = array(0, dim = numOfMCMC)
  dpg = array(0, dim = numOfMCMC)
  dd = array(0, dim = c(dims, dims, numOfMCMC))

  if (is.null(tssb)) {
    root <- Normal$new(...)
    tssb <- TssbMCMC$new(root, data = data)
  } else {
    root <- tssb$root$node
  }

  for (i in burnIn:numOfMCMC) {

    tssb$ResampleAssignments()

    tssb$CullTree()

    tssb$ResampleNodeParameters()

    root$ResampleHyperParams()

    tssb$ResampleSticks()

    tssb$ResampleStickOrders()

    tssb$ResampleHypers()

    res <- tssb$GetLogMarginalDataLikelihood()

     if (i >0) {
      ll[i] <- res$ll
      pp[,,i] <- Reduce(rbind, Map(function(x) x$params, tssb$assignments))
      ss[,,i] <- Reduce(rbind, Map(function(x) x$sigma, tssb$assignments))
      dpa[i] <- tssb$dpAlpha
      dpl[i] <- tssb$dpLambda
      dpg[i] <- tssb$dpGamma
      dd[,,i] <- root$GetDrift()
      gg[[i]] <- tssb$ConvertTssbToIgraph()$g
    }

    if (i%%100 == 0) {
      cat("iteration:", i,  "average datum marginal likelihood:", res$ll/numOfData,
          "nodes:", length(res$ww), "big nodes:",  length(res$ww > 0.01),
          "\n")
    }
  }

  return(list(tssb = tssb,
              likelihood = ll,
              mu = pp,
              sigma = ss,
              dpAlpha = dpa,
              dpLambda= dpl,
              dpGamma = dpg,
              igraph = gg,
              drift = dd))
}