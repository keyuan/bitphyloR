n0 <- Node$new()
tssb <- TSSB$new(n0, data = matrix(rnorm(50),50,1), dpAlpha = 1, dpGamma = 1/5,
                 dpLambda = 1/2)
res <- tssb$GetMixture()

res1 <- tssb$ConvertTssbToIgraph()

TestDescend <- function(root){
  Descend <- function(root) {
    if (!is.null(root$sticks)) {
      root$main <- length(root$children)
      for (i in 1:length(root$children)) {
        Descend(root <<- root$children[[i]])
      }
    }
    return(root)
  }
  return(Descend(root<<-root))
}
TestDescend(root<<-tssb$root)