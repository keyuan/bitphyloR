library(microbenchmark)

#set.seed(9)
# n0 <- Node$new()
# tssb <- TSSB$new(n0, data = matrix(rnorm(50),50,1), dpAlpha = 1, dpGamma = 1,
#                  dpLambda = 1)
# res1 <- tssb$GetMixture()
# res2 <- tssb$ConvertTssbToIgraph()
#
# plot(res2$g, layout = layout.reingold.tilford(res2$g),
#      vertex.label = get.vertex.attribute(res2$g, name = "size"),
#      vertex.size = 15, edge.arrow.size = 0.5)
#
# tt <- tssb$root
#
# tssb$CullTree()$CullTree()
# tt1 <- tssb$root
# res3 <- tssb$ConvertTssbToIgraph()
# plot(res3$g,
#      layout = layout.reingold.tilford(res3$g),
#      vertex.label = get.vertex.attribute(res3$g, name = "size"),
#      vertex.size = 15, edge.arrow.size = 0.5)

set.seed(10)
n = 1
res1 = vector("numeric", n)
res2 = res1
res3 = res1
res4 = res1
res5 = res1
res6 = res1
for (i in seq_along(res1)){
  cat(i)
  n0 <- Node$new()
  tssbMCMC <- TssbMCMC$new(n0,
                           data = matrix(rnorm(50),50,1),
                           dpAlpha = 1,
                           dpGamma = 1,
                           dpLambda = 1)
  ww <- tssbMCMC$GetMixture()$weight
  res1[i] <- sum(ww)
  res2[i] <- length(ww)
#   root1 <- tssbMCMC$root
#   g <- tssbMCMC$ConvertTssbToIgraph()$g
#   plot(g,
#        layout = layout.reingold.tilford(g),
#        vertex.label = get.vertex.attribute(g, name = "size"),
#        vertex.size = 15, edge.arrow.size = 0.5)

  tssbMCMC$CullTree()
  ww <- tssbMCMC$GetMixture()$weight
  res3[i] <- sum(ww)
  res4[i] <- length(ww)
#   root2 <- tssbMCMC$root
#   g <- tssbMCMC$ConvertTssbToIgraph()$g
#   plot(g,
#        layout = layout.reingold.tilford(g),
#        vertex.label = get.vertex.attribute(g, name = "size"),
#        vertex.size = 15, edge.arrow.size = 0.5)

  ww <- tssbMCMC$ResampleSticks()$GetMixture()$weight
  res5[i] <- sum(ww)
  res6[i] <- length(ww)

}





