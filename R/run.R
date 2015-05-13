#' Fit a tssb with normal node with data via mcmc
#' @param data a data matrix
#' @param numOfMCMC mcmc iterations
#' @return A list of traces
#' @export

RunNormal <- function (data, numOfMCMC, tssb = NULL) {

  priorSigmaScale <- var(data)
  priorDriftScale <- priorSigmaScale
  rootNode <- Normal(priorSigmaScale = priorSigmaScale, priorDriftScale = priorDriftScale)
  tssb <- TssbMCMC$new(rootNode, data = data)
  dims <- ncol(data)
  numOfData <- nrow(data)

  gg = vector(mode = "list", length = numOfMCMC)
  ll = array(0, dim = numOfMCMC)
  pp = array(0, dim = c(numOfMCMC, numOfData))
  ss = array(0, dim = c(numOfMCMC, numOfData))
  dpa = array(0, dim = numOfMCMC)
  dpl = array(0, dim = numOfMCMC)
  dpg = array(0, dim = numOfMCMC)
  dd = array(0, dim = numOfMCMC)


  for (n in seq_len(numOfMCMC)) {
    tssbMCMC$ResampleAssignments()

    tssbMCMC$CullTree()

    tssbMCMC$ResampleNodeParameters()

    rootNode$ResampleHyperParams()

    tssbMCMC$ResampleSticks()

    tssbMCMC$ResampleStickOrders()

    tssbMCMC$ResampleHypers()

    ll[i] <- tssbMCMC$GetLogMarginalDataLikelihood()
    pp[i, ] <- sapply(tssbMCMC$assignments, function(x) x$params)
    ss[i, ] <- sapply(tssbMCMC$assignments, function(x) x$sigma)

    dpa[i] <- tssbMCMC$dpAlpha
    dpl[i] <- tssbMCMC$dpLambda
    dpg[i] <- tssbMCMC$dpGamma
    dd[i] <- n0$GetDrift()
    gg[[i]] <- tssbMCMC$ConvertTssbToIgraph()$g

  }

  return(list(likelihood = ll,
              mu = pp,
              sigma = ss,
              dpAlpha = dpa,
              dpLambda= dpl,
              dpGamma = dpg,
              igraph = gg))
}