# metrics


brier <- function(x, y){
  N <- length(x)
  bs <- (1/N) * sum((x-y)^2)
  return(bs)
}

brier_cal <- function(x, y){
  f <- unique(x)
  k <- length(f)
  N <- length(x)

  tot <- 0
  for(i in 1:k){
    inds <- which(x == f[i])
    nk <- length(x[inds])
    okb <- sum(y[inds])/nk
    tot <- tot + nk * (f[i] - okb)^2
  }

  return((1/N) * tot)
}

brier_ref <- function(x, y){
  f <- unique(x)
  k <- length(f)
  N <- length(x)

  tot <- 0
  for(i in 1:k){
    inds <- which(x == f[i])
    nk <- length(x[inds])
    okb <- sum(y[inds])/nk
    tot <- tot + nk * okb * (1 - okb)  # we will always have zero refinement if we don't have repeat x's....
  }

  return((1/N) * tot)
}


