library(igraph)
library(RBGL) ## masks functions degree() edges() from graph and transitivity() from igraph
## contains edmondsOptimumBranching()
library(Matrix)
rm(list = ls())

## choose root
setRoot <- function(W){
    sumIn <- apply(W, 2, sum)
    root <- which.min(sumIn)
    W[,root] <- 0
    diag(W) <- 0
    return(list(W=W, root=root))
}

setMinInEdgeGraph <- function(W){
    n <- dim(W)[1]
    diag(W) <- -1/0
    selectedInEdges <- max.col(t(W))
    newW <- matrix(0,n,n)
    for (i in 1:n){
        newW[selectedInEdges[i],i] <- W[selectedInEdges[i],i]
    }
    return(newW)
}

## Scales scores so that all entries of the matrix are >= 0 , where 0 indicates no edge
scaleScore <- function(score) {
  score <- score - min(score) + 1
  diag(score) <- 0
  return(score)
}

## Finds optimum branching, given a weighted adjacency matrix. (Optimum = maximum)
EdmondsOpt <- function(score,ret='graph') {
  if (all(score==t(score))) {
    stop("only appropriate for directed graphs")
  }
  if (all(score<=0)) {
    score <- scaleScore(score)
  }
  nv <- nrow(score) ## number of vertices
  ne <- prod(dim(score)) - nv
  em <- which(!diag(nrow(score)),arr.ind=TRUE)
  em <- rbind(from=em[order(em[,1]),1],to=em[order(em[,1]),2])
  eW <- score[t(em)]
  ans <- .Call("edmondsOptimumBranching", as.integer(nv), as.integer(ne),
               as.integer(em - 1), as.double(eW), PACKAGE = "RBGL")
  el <- t(ans[[1]]+1) ## start node numbering with 1
  eW <- ans[[2]]
  if (!is.null(rownames(score))) {
    el <- matrix(rownames(score)[el],ncol=2,byrow=FALSE)
  }
  if (ret=='graph'){
    el <- el[order(el[,1]),] ## order nodes to make sure plotting function works
    g <- graph.edgelist(el,directed = TRUE)
    return(g)
  } else if (ret=='adjMat') {
    adjMat <- matrix(0,nrow=nv,ncol=nv)
    rownames(adjMat) <- colnames(adjMat) <- rownames(score)
    adjMat[el] <- 1
    return(adjMat)
  }
}

plotG <- function(W){
    plot(graph.adjacency(W, weighted = T))
}

plotT <- function(W){
    root <- which(colSums(W) == 0)
    g <- graph.adjacency(W, weighted=T)
    l <- layout.reingold.tilford(g,root=root)
    plot(g, layout=l)
}

boundProb <- function(x){
    return( (1.0-.Machine$double.eps)*(x-0.5)+0.5 )
}


getWeights <- function(parentStates, ratemat, branchLength){
    transProb <- expm(ratemat * branchLength)
    m1 <- transProb[2,2]
    if (m1 == 1.0){
        m1 <- boundProb(m1)
    }
    m2 <- transProb[1,2]
    return(parentStates*m1 + (!parentStates)*m2)
}


getLocalLikhd <- function(childSeq, parentSeq, ratemat, t){
    prob <- getWeights(parentSeq, ratemat, t)
    prob <- childSeq*prob + (!childSeq)*(1-prob)
    return( sum(log(prob)) )
}

getMaxLocalLikhd <-function(childSeq, parentSeq, ratemat, tspace=NULL){
    if (is.null(tspace)){
        tspace <- seq(0,1,length.out=100)
    }
    ll <- sapply(tspace, function (i) getLocalLikhd(childSeq,
        parentSeq, ratemat, i) )
    optLL <- max(ll)
    optT <- tspace[ which.max(ll) ]
    return( list(optLL = optLL, optT = optT) )
}

getPhyloDist <- function(seq, ratemat){
    n <- dim(seq)[1]
    likhdW <- matrix(0,n,n)
    branchW <- matrix(0,n,n)
    for (i in 1:n){
        for (j in 1:n){
            if (i !=j){
                maxLocalLikhd <- getMaxLocalLikhd(seq[i,], seq[j,], ratemat)
                likhdW[j,i] <- maxLocalLikhd$optLL
                branchW[j,i] <- maxLocalLikhd$optT
            }
        }
    }
    return(list(likhdW = likhdW, branchW = branchW))
}

getPhyloDist2 <- function(seq, ratemat){
  n <- dim(seq)[1]
  likhdW <- matrix(0,n,n)
  branchW <- matrix(0,n,n)
  loopMat <- combn(n, 2)
  loopMat <- t(cbind(loopMat, loopMat[c(2,1),]))
  res <- foreach(i=1:(n*(n-1)),
                 .packages = 'Matrix',
                 .export='getMaxLocalLikhd') %do% {
    index <- loopMat[i,]
    getMaxLocalLikhd(seq[index[2],], seq[index[1],], ratemat)
  }
  res <- unlist(res)
  likhdW[loopMat] <- res[which(names(res) == 'optLL')]
  branchW[loopMat] <- res[which(names(res) == 'optT')]

  return(list(likhdW = likhdW, branchW = branchW))
}

getRateMatrix <- function(data, mode='methyl'){
    mrate <- mean(rowSums(data, na.rm = TRUE) / rowSums(!is.na(data)))
    urate <- 1-mrate
    ratemat <- matrix(0, 2, 2)
    ratemat[1,] <- c(-mrate, mrate)
    if (mode == 'methyl'){
        scale <- 2/(3*urate*mrate)
        ratemat[2,] <- c(urate/2, -urate/2)
        ratemat <- scale*ratemat
    }else{
        scale <- 1/(urate*mrate)
        ratemat <- scale*ratemat
    }
    return(ratemat)
}

getDMST <- function(data, mode='methyl'){
    ratemat <- getRateMatrix(data, mode)
    distmats <- getPhyloDist(data, ratemat)
    likhdW <- distmats$likhdW
    branchW <- distmats$branchW
    A <- EdmondsOpt(likhdW, ret ='adjMat')
    mstW1 <- distmats$likhdW*A
    mstW2 <- distmats$branchW*A
    return(list(W1=mstW1, W2 = mstW2,  distmats=distmats, ratemat=ratemat))
}



f1 <- '~/Dropbox/cancer-evolution/tssb/data/full_methy/noisy_full_methy_8_2000_genotype.csv'
f2 <- '~/Dropbox/cancer-evolution/tssb/data/noisy-full-obs/noisy_fullsyn_8_2000_genotype.csv'

data1 <- t(read.csv(f1, as.is=TRUE))
data2 <- t(read.csv(f2, as.is=TRUE))

tic <- proc.time()
mst1 <- getDMST(data1, mode = 'methyl')
mst2 <- getDMST(data2, mode = 'mut')
toc <- proc.time()-tic

stopCluster(cl)
