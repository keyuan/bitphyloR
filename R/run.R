RunNormal <- function () {

  priorSigmaScale <- var(data)
  priorDriftScale <- priorSigmaScale
  rootNode <- Normal(priorSigmaScale = priorSigmaScale, priorDriftScale = priorDriftScale)
  tssb <- TssbMCMC$new(rootNode, data = data)




  for (n in seq_len(numOfMCMC)) {
    tssbMCMC$ResampleAssignments()

    tssbMCMC$CullTree()

    tssbMCMC$ResampleNodeParameters()

    rootNode$ResampleHyperParams()

    tssbMCMC$ResampleSticks()

    tssbMCMC$ResampleStickOrders()

    tssbMCMC$ResampleHypers()
  }





}