llike <- function(x, n, delta2, Delta){
  return(-1/2 * (n * log(2 * pi) + sum(log(Delta + x)) + sum(delta2 / (Delta +  x))))
}

dllike <- function(x, delta2, Delta){
  return(-0.5 * (sum(1 / (Delta + x)) - sum(delta2 / (Delta + x)^2)))
}
d2llike <- function(x, delta2, Delta){
  return(-0.5 * (-sum(1 / (Delta + x)^2) + 2 * sum(delta2 / (Delta + x)^3)))
}

sqrt_ginv <- function(X, tol = sqrt(1e-16)){
  Xsvd <- svd(X);
  u <- Xsvd$u; d <- Xsvd$d; v <- t(Xsvd$v)    
  pos <- d > max(tol * d[1], 0)
  res <- t(v)[,pos] %*% sqrt(diag(1 / d[pos])) %*% t(u[, pos])
  return(res)
}

newton_optim <- function(x, delta2, Delta, maxit = 1e4, tol = 1e-6){
  i <- 0
  while(abs(dllike(x, delta2, Delta)) > tol & i < maxit){
    x <- x - dllike(x, delta2, Delta) / d2llike(x, delta2, Delta)
    i <- i + 1
  }
  return(x)
}
