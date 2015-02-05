n0 <- Node$new()
tssb <- TSSB$new(n0, data = matrix(rnorm(50),50,1), dpAlpha = 1, dpGamma = 1,
                 dpLambda = 1)
res <- tssb$GetMixture()

expect_true((1-sum(res$weight)) < 0.3)
