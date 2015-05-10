library(microbenchmark)

set.seed(10)
n = 1
res1 = vector("numeric", n)
res2 = res1
res3 = res1
res4 = res1
res5 = res1
res6 = res1

n0 <- Normal$new()
tssbMCMC <- TssbMCMC$new(n0,
                         data = matrix(rnorm(50),50,1),
                         dpAlpha = 1,
                         dpGamma = 1,
                         dpLambda = 1)


tssbMCMC$ResampleAssignments()

tssbMCMC$CullTree()

tssbMCMC$ResampleNodeParameters()

tssbMCMC$ResampleSticks()

tssbMCMC$ResampleStickOrders()

tssbMCMC$CullTree()

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

fx <- function(x) {
  - 0.5*(x-200)^2/(10)
}

n = 1000
x = array(0, dim = n)
for (i in seq_len(n-1)) {
  x[i+1] = SliceSampler(x[i], fx, compwise = F, stepOut = T)
}
