n0 <- Node$new()
tssb <- TSSB$new(n0, data = matrix(rnorm(50),50,1), dpAlpha = 1, dpGamma = 1/5,
                 dpLambda = 1/2)
res <- tssb$GetMixture()

res1 <- tssb$ConvertTssbToIgraph()

plot(res1$g)