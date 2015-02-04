

tssb <- TSSB$new(n0, data = matrix(rnorm(50),10000,1))
res <- tssb$GetMixture()
(1-sum(res$weights) < 0.2)