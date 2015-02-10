n0 <- Node$new()
tssb <- TSSB$new(n0, data = matrix(rnorm(50),50,1), dpAlpha = 1, dpGamma = 1/5,
                 dpLambda = 1/2)
res1 <- tssb$GetMixture()

sapply(res1$node, function(x) x$GetNumOfLocalData())

res2 <- tssb$ConvertTssbToIgraph()

tt <- tssb$root


Descend1 <- function(root) {
  res <- unlist(Map(Descend1, root$children), recursive = F, use.names = F)
  if (length(res) == 0) {
    return(list(counts = root$node$GetNumOfLocalData(), root = root))
  }
  counts <- unlist(res[seq(1, length(res), by = 2)])
  root$children <- res[seq(2, length(res), by = 2)]
  keep <- length(TrimZeros(counts, trim = "b"))

  if (keep == 0) {
    root$sticks <- NULL
    root$children <- list()
  } else {
    root$sticks <- root$sticks[1:keep]
    root$children <- root$children[1:keep]
  }

  return(list(counts = sum(counts) + root$node$GetNumOfLocalData(),
              root = root))
}

td <- Descend1(tssb$root)

