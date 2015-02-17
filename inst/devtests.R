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

res1 = c()
res2 = c()
res3 = c()
res4 = c()
res5 = c()
res6 = c()
for (i in 1:1){
  print(i)
  n0 <- Node$new()
  tssbMCMC <- TssbMCMC$new(n0,
                           data = matrix(rnorm(50),50,1),
                           dpAlpha = 1,
                           dpGamma = 1,
                           dpLambda = 1)
  ww <- tssbMCMC$GetMixture()$weight
  res1 <- c(res1, sum(ww))
  res2 <- c(res2, length(ww))
  root1 <- tssbMCMC$root

  tssbMCMC$CullTree()
  ww <- tssbMCMC$GetMixture()$weight
  res3 <- c(res3, sum(ww))
  res4 <- c(res4, length(ww))
  root2 <- tssbMCMC$root

  ww <- tssbMCMC$ResampleSticks()$GetMixture()$weight
  res5 <- c(res5, sum(ww))
  res6 <- c(res6, length(ww))
  root3 <- tssbMCMC$root

}






