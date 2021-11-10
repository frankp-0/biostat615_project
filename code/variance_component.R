run_vc_optimizer <- function(x, sUt, Ce, ind){
    b <- t(sUt) %*% x[, ..ind]
    ind1 <- ind + 1
    se <- x[, ..ind1]
    K = t(sUt) %*% diag(se) %*% Ce %*% diag(se) %*% sUt
    return(vcm_optimization(b, K))
}

sqrt_ginv <- function(X, tol = sqrt(2.22044604925e-16)){
    Xsvd <- svd(X);
    u <- Xsvd$u; s <- Xsvd$d; vh <- Xsvd$v
    pos <- s > max(tol * s[1], 0)
    if (all(pos)){
        res <- t(vh) %*% sqrt(diag(1 / s)) %*% t(u)
    } else if (all(!pos)){
        res <- matrix(0, nrow = nrow(X), ncol = ncol(X))
    } else {
        res <- t(vh)[, pos] %*% sqrt(diag(1 / s[pos])) %*% t(u[, pos])
    }
    return(res)
}

vcm_optimization <- function(){
    }

LL_fun <- function(x, n, sq, w){
    return(-0.5 * (n * log(2 * pi) + sum(log(w + x)) + sum(sq / (w +  x))))
}

LLp_fun <- function(x, sq, w){
    return(-0.5 * (sum(1 / (w + x)) - sum(sq / (w + x)^2)))    
    }

LLdp_fun <- function(x, sq, w){
    return(-0.5*(-sum(1/(w+x)**2)+2*sum(sq/(w+x)**3)))    
    }

NR_root <- function(f, df, x, sq, w, i = 0, iter_max = 10000, tol = sqrt(2.22044604925e-16)){
    while(abs(f(x, sq, w)) > tol){
        x <- x - f(x, sq, w) / df(x, sq, w)
        i <- i + 1
        if (i == iter_max){break}
    }
    return(x)
}

vcm_optimization <- function(b, K, tol = sqrt(2.22044604925e-16)){
    n = length(b)
    eigk <- eigen(K)
    l <- eigk$values[n:1]; cv <- -eigk$vectors[,n:1]
    p = l > tol 
    if(all(p)){
        crossP = t(cv) %*% (b)
        sq = crossP**2;
        w = l;
    } else{
        cr = t(cv[,p]) %*% b;
        sq = cr**2;
        w = l[p]
    }
    t <- 10**((-36:23)/4)
    init <- it[which.max(sapply(it, function(i) LL_fun(i, n, sq, w)))]
    tausq <- NR_root(LLp_fun, LLdp_fun, init, sq, w)
    if (tausq <0){tausq = 0}
    null_ll = LL_fun(0, n, sq, w)
    alt_ll = LL_fun(tausq, n, sq, w)
    if(alt_ll < null_ll){tausq = 0; alt_ll = null_ll}
    pleio_stat <- -2 * (null_ll - alt_ll)
    return(list(tausq = tausq, pleio_stat = pleio_stat))
    }
l    
