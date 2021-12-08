source("helpers.R")

vc_test <- function(X, Omega, C_e, tol = 1e-8){
  omega_inv_sq <- sqrt_ginv(Omega)
  n <- length(X) / 2
  ind <- seq(1, 2 * n, 2)
  eta_prime <- t(omega_inv_sq) %*% X[ind]
  se <- X[ind + 1]
  eig <- eigen(t(omega_inv_sq) %*% diag(se) %*% C_e %*% diag(se) %*% omega_inv_sq,
               symmetric = T)
  eigval <- eig$values[n:1]; eigvec <- -eig$vectors[,n:1]
  delta2 <- (t(eigvec[, eigval > tol]) %*% eta_prime)^2;
  Delta <- eigval[eigval > tol]
  it <- 10^((-36:23)/4)
  init <- it[which.max(sapply(it, function(i) llike(i, n, delta2, Delta)))]
  tausq <- max(0, newton_optim(init, delta2, Delta))
  ll0 <- llike(0, n, delta2, Delta)
  ll1 <- llike(tausq, n, delta2, Delta)
  if(ll1 < ll0){
    tausq <- 0
    ll1 <- ll0
  }
  stat <- -2 * (ll0 - ll1)
  return(c(tausq, stat))
}
