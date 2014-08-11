## EM for Bernoulli and logistic mixtures
## Author: Ke Yuan
library(Matrix)

q = matrix(0,2,2)
q[1,] = c(-1,1)
tt = c(0.1, 0.2, 0.3, 0.3, 0.5, 0.5)
P = lapply(1:6, function(i) expm(q*tt[i]) )
n = 20
root = matrix(0,1,n)

desend <- function(parent, P, n){
  weights <- parent * P[2,2] + (!parent) * P[1,2]
  return(runif(n) < weights)
}

logSumExp <- function(x, dim=1){
  ## todo: Bugs when dim = 2
  x <- as.matrix(x)
  nargin <- length(as.list(match.call())) - 1
  if (nargin == 1){
      rdim <- dim(x)[2:1]
      dim <- which(rdim!=1)[1]
  }
  y <- apply(x, dim, max)
  x <- x - y
  s <- y + log(apply(exp(x), dim, sum))
  return(s)
}

repmat <- function(x, n, m){
  x <- as.matrix(x)
  return(kronecker(matrix(1,n,m),x))
}

sigmoid <- function(x){
  return(1/(1+exp(-x)))
}

sigmoidln <- function(x){
  return(-log(1+exp(-x)))
}

constriant <- function(x){
  if (length(which(x > 0.9999) ) >0 ){
    x[which(x>0.9999)] = 0.9999
  }
  if (length(which(x < 1e-20)) > 0){
    x[which(x<1e-20)] = 1e-20
  }
  return(x)
}

logit <- function(x){
  x <- constriant(x)
  return(-log(1/x-1))
}

logit_llh <- function(x, y){
  x = t(as.matrix(x))
  x = repmat(x, dim(y)[1], 1)
  llh = y*sigmoidln(x) + (1-y)*sigmoidln(-x)
  return(rowSums(llh))
}

bernoulli_llh <- function(x,y){
  x = t(as.matrix(x))

  x = repmat(x, dim(y)[1], 1)
  llh = y*log(x) + (1-y)*log(1-x)
  return(rowSums(llh))
}

## EM

ExpectMax <- function(data, K, func, eps, totem){

  ## Initial
  P <- dim(data)[2]
  N <- dim(data)[1]
  Eparam <- matrix(runif(K*P), K, P)
  a <- rep(0.5, K)
  Epi <- rdirichlet(1, a)
  vecEpi <- repmat(Epi, N, 1)
  veclogEpi <- log(vecEpi)
  fxln <- matrix(0, N, K)
  llh <- rep(0,totem)

  for (nem in 1:totem){

    ## E-step: compute responsabilies
    for (k in 1:K){
      fxln[,k] <- do.call(func, list(Eparam[k,], data))
    }
    llh_tmp <- logSumExp(veclogEpi + fxln, 1)
    llh_tmp <- repmat(llh_tmp, 1, K)
    resln <- veclogEpi + fxln - llh_tmp

    ## E-step: compute log-marginal
    llh[nem] <- sum(llh_tmp[,1])
    res <- exp(resln)


    ## M-step: compute node params
    Nk <- colSums(res)
    vecNk <- repmat(Nk, 1, P)
    data <- as.matrix(data)
    xkbar <- t(res)%*%data/vecNk

    if (func == 'logit_llh'){
      Eparam <- logit(xkbar)
    }
    else {
      Eparam <- constriant(xkbar)
    }

    ## M-step: compute weights
    Epi <- Nk/N
    vecEpi <- t(repmat(Epi, 1, N))
    veclogEpi <- log(vecEpi)

    ## stopping condition
    if (nem>2){
      gain = llh[nem] - llh[nem-1]
      cat(sprintf("\rEM steps: %i, log-marginal likelihood: %f, relative gain: %f\r",
                    nem, llh[nem], abs(gain/llh[nem-1])))
      if (abs( gain/llh[nem-1] ) < eps && gain > 0) {
        cat(sprintf('\nconverged.\n'))
        break
      }
    }
  }

  pred_label <- apply(res, 1, which.max)

  return( list(pred_label=pred_label, llh=llh[1:nem], Eparam=Eparam, Epi=Epi) )

}


filepath <- '' # add your own path
files <- dir(filepath, pattern = 'mutmat')

for (file in files){

  cat(sprintf('Processing %s ...\n', file))
  data <- read.csv(paste(filepath,file,sep=''))

  tmpname <- substr(file, 22, nchar(file)-11)

  true_label <- data[,9]
  data <- data[,-9]
  K <- 7
  eps <- 1e-6
  totem <- 1000

  cat(sprintf('Fitting Bernoulli mixture:\n'))
  em_bernoulli <- ExpectMax(data, K, 'bernoulli_llh', eps, totem)
  cat(sprintf('Fitting logistic mixture:\n'))
  em_logit <- ExpectMax(data, K, 'logit_llh', eps, totem)

  pred_label_bernoulli <- em_bernoulli$pred_label
  pred_label_logit     <- em_logit$pred_label

  labels <- data.frame( bernoulli = pred_label_bernoulli,
                        logit = pred_label_logit,
                        truth = true_label)

  write.csv( labels, file = paste('simple_mixtures_labels_', tmpname, '.csv',
                                  sep=''), row.names=F )
}


