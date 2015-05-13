library(microbenchmark)

set.seed(10)

m = 100
testData = rbind(matrix(0.051*rnorm(m)+0.4,m,1), matrix( 0.01*rnorm(m)+0.6,m,1))
n0 <- Normal$new(priorSigmaScale = diag(0.1, length(0.1)))
tssbMCMC <- TssbMCMC$new(n0,
                         data = testData,
                         dpAlpha = 1,
                         dpGamma = 1,
                         dpLambda = 1)

#
n = 1000
gg = vector(mode = "list", length = n)
ll = array(0, dim = n)
pp = array(0, dim = c(n, 2*m))
ss = array(0, dim = c(n, 2*m))
dpa = array(0, dim = n)
dpl = array(0, dim = n)
dpg = array(0, dim = n)
dd = array(0, dim = n)

for (i in seq_len(n)) {

  if (i%%50== 0) {
    print(i)
  }

  tssbMCMC$ResampleAssignments()

  tssbMCMC$CullTree()

  tssbMCMC$ResampleNodeParameters()

  n0$ResampleHyperParams()

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


g <- tssbMCMC$ConvertTssbToIgraph()$g

for (i in 1:n) {
  plot(gg[[i]],
       layout = layout.reingold.tilford(gg[[i]]),
       vertex.label = get.vertex.attribute(gg[[i]], name = "size"),
       vertex.size = 15, edge.arrow.size = 0.5)

}

## tssbMCMC$ResampleNodeParameters()
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





