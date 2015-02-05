n0 <- Node$new()
tssb <- TSSB$new(n0, data = matrix(rnorm(50),50,1), dpAlpha = 1, dpGamma = 1,
                 dpLambda = 1)
res <- tssb$GetMixture()

res1 <- tssb$ConvertTssbToIgraph()
