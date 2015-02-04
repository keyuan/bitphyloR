SticksToEdges <- function(sticks){
  return(1.0 - cumprod(1.0 - sticks))
}

