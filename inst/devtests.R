library(microbenchmark)

set.seed(10)

testData = matrix(rnorm(50),50,1)
n0 <- Normal$new()
tssbMCMC <- TssbMCMC$new(n0,
                         data = testData,
                         dpAlpha = 1,
                         dpGamma = 1,
                         dpLambda = 1)


n = 100

ll = array(0, dim = n)

for (i in seq_len(n)) {

  print(i)

  tssbMCMC$ResampleAssignments()


  tssbMCMC$CullTree()

  tssbMCMC$ResampleNodeParameters()


  tssbMCMC$ResampleSticks()


  tssbMCMC$ResampleStickOrders()


  tssbMCMC$CullTree()


  #tssbMCMC$ResampleHypers()

  ll[i] <- tssbMCMC$GetLogMarginalDataLikelihood()

}



#   root1 <- tssbMCMC$root
#   g <- tssbMCMC$ConvertTssbToIgraph()$g
#   plot(g,
#        layout = layout.reingold.tilford(g),
#        vertex.label = get.vertex.attribute(g, name = "size"),
#        vertex.size = 15, edge.arrow.size = 0.5)

#   root2 <- tssbMCMC$root
#   g <- tssbMCMC$ConvertTssbToIgraph()$g
#   plot(g,
#        layout = layout.reingold.tilford(g),
#        vertex.label = get.vertex.attribute(g, name = "size"),
#        vertex.size = 15, edge.arrow.size = 0.5)

# fx <- function(x, y) {
#   print(y)
#   - 0.5*(x-0)^2/(1) - 0.5*(y-0)^2/(1)
# }
#
# n = 1000
# x = array(0, dim = c(n,1))
# for (i in seq_len(n-1)) {
#   x[i+1,] = SliceSampler(x[i,], fx, compwise = T, stepOut = T, y = 10)
# }





