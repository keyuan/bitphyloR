
ConvertCellPrevToAllelePrev <- function (mode = "pyclone", ... ) {
  if (mode == "pyclone") {
    ap = MapPyClone(...)
  }
}

MapPyClone <- function(x,t,cn,cr,cv,epi,bv){
  p1 <- (1-t)*cn
  p2 <- t*(1-x)*cr
  p3 <- t*x*cv
  p <- p1+p2+p3
  u1 <- epi
  u2 <- epi
  if (bv == cv) {
    u3 <- 1-epi
  } else {
    u3 <- bv/cv
  }
  y <- p1*u1+p2*u2+p3*u3
  return(y/p)
}


PyCloneLikelihood <- function() {

}