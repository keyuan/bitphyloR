library(microbenchmark)

set.seed(10)
#
m <- 200
dims <- 2
testData <- rbind(matrix(0.051*rnorm(m)+0.4,m/2,dims),
                  matrix( 0.01*rnorm(m)+0.6,m/2,dims))

empCov <- cov(testData)
priorSigmaScale = empCov + diag(rep(1e-6,dims), dims)
priorDriftScale = priorSigmaScale

traces <- RunNormal(data = testData,
                    priorSigmaScale = priorSigmaScale,
                    priorDriftScale = priorDriftScale)


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





