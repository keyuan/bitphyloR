test_that("classname contains 'TSSB'", {
  n0 <- Node$new()
  tssb <- TSSB$new(n0, data = matrix(rnorm(50),50,1), dpAlpha = 1, dpGamma = 1,
                   dpLambda = 1)
  expect_that("TSSB" %in% class(tssb), is_true())
})

test_that("Sum of sticks needs to be close to one", {
  n0 <- Node$new()
  tssb <- TSSB$new(n0, data = matrix(rnorm(50),50,1), dpAlpha = 1, dpGamma = 1,
                   dpLambda = 1)
  res <- tssb$GetMixture()
  expect_true((1-sum(res$weight)) < 0.3)
})




