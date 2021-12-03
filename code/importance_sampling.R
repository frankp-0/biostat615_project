N <- 100
gwas_N <- read.csv("~/pleio/N_gwas.csv")[,-1]
U <- as.matrix(read.csv("~/pleio/Sg.csv")[,-1])
Ce <- as.matrix(read.csv("~/pleio/Ce.csv")[-1])
importance_sampling <- function(N, gwas_N, U, Ce, tol = 2.22044604925e-16**0.5, out){
    se <- 1/sqrt(gwas_N)
    nstudy <- n <- length(se)
    D <- diag(se) %*% Ce %*% diag(se)
    null_D <- diag(rep(1, n)) %*% Ce %*% diag(rep(1, n))
    Uinv_sqrt <- sqrt_ginv(U)
    K <- t(Uinv_sqrt %*% D %*% Uinv_sqrt)
    eigk <- eigen(K)
    w <- rev(eigk$values)
    v <- eigk$vectors[,18:1]
    t_v <- t(v)
    pos <- w > tol
    w_pos <- w[pos]
    t_v_pos <- t_v[,pos]
    c_Pj <- c(1,1.1,1.2,1.3,1.4,1.7,2,2.5,3,4,5)
    nPj = length(c_Pj); mean_P = rep(0, nPj); alpha = rep(1/nPj, nPj)
    cov_c <- lapply(1:nPj, function(i) {
        diag(rep(c_Pj[i], n)) %*% null_D %*% diag(rep(c_Pj[i], n))
    })
    mean_P <- lapply(1:nPj, function(i) {
        rep(mean_P[i], n)
    })
    P <- list(means = mean_P, cov = cov_c)
    Ppdf <- mixture_sampling(N, alpha, P)
    eta_df <- sweep(Ppdf, 2, se, `*`)
    transformed_df <- eta_df %*% t(Uinv_sqrt)
    data <- apply(transformed_df, 1, function(b) {
        vcm_optimization_IS(b, n, w, t_v)
    })
    Sdelpy = data
    d_Q <- mvtnorm::dmvnorm(Ppdf, rep(0, n), Ce)
    d_P <- lapply(1:length(P$means), function(i) {
        mvtnorm::dmvnorm(Ppdf, P$means[[i]], P$cov[[i]])
    })
    thres_vec <- c(0, 40^(seq(-5,0,length=20))[1:19], 40^(seq(0,1,length=20)))
    Palpha <- rowSums(t(do.call("rbind", d_P)) * alpha)
    pvalues <- sapply(thres_vec, function(thres) thres_estimate_pvalue(thres, Sdelpy, Palpha, alpha, d_Q, d_P, nPj, N))
    pvalue_df <- data.frame(thres = thres_vec, p = pvalues)
    pvalue_df <- pvalue_df[order(pvalue_df$p, decreasing = T),]
    write.csv(pvalue_df, file = out)
}

thres_estimate_pvalue <- function(thres, Sdelpy, Palpha, alpha, d_Q, d_P, nPj, N){
    h <- as.numeric(Sdelpy > thres)
    m <- h * d_Q / Palpha
    cov_tm <- sapply(1:length(alpha), function(i) cov(cbind(m, d_P[[i]]))[1,2])
    cov_t <- cov(t(do.call("rbind", d_P))/Palpha)
    inv_cov_t <- svd_inv(cov_t)
    denominator <- rowSums(t(do.call("rbind", d_P) * alpha))
    betas <- (inv_cov_t %*% cov_tm)[,1]
    control_variate <- rowSums(t(do.call("rbind", d_P) * betas))
    nominator <- h * d_Q - control_variate
    IS_estim <- sum(nominator / denominator) / N + sum(betas)
    return(IS_estim)
}

svd_inv <- function(cov_t){
    s <- svd(cov_t)
    ds <- diag(1/s$d[1:(length(s$d)-1)])
    us <- s$u[,-ncol(s$u)]
    vs <- s$v[,-ncol(s$v)]
    return(vs %*% ds %*% t(us))
    }

mixture_sampling <- function(N, alpha, P){
    K <- length(alpha)
    counts <- rep(0, K)
    comp <- sample(1:K, N, replace = T, prob = alpha)
    counts[as.integer(names(table(comp)))] <- table(comp)
    P.pdf <- lapply(which(counts != 0), function(i) {
        MASS::mvrnorm(n = counts[i], mu = rep(P$means[[i]], ), Sigma = P$cov[[i]])
    })
    return(do.call("rbind", P.pdf))
}

vcm_optimization_IS <- function(b, n, w, t_v, tol = 1e-17){
    tt <- 10^((-36:24)/4)
    crossP <- t_v %*% b
    P_sq <- crossP^2
    init <- tt[which.max(sapply(tt, function(i) LL_fun(i, n, P_sq, w)))]
    mle_tausq <- NR_root(LLp_fun, LLdp_fun, init, P_sq, w)
    mle_tausq <- max(0, mle_tausq)
    null_ll <- LL_fun(0, n, P_sq, w)
    alt_ll <- LL_fun(mle_tausq, n, P_sq, w)
    if(alt_ll < null_ll){
        mle_tausq <- 0
        alt_ll <- null_ll
    }
    return(-2 * (null_ll - alt_ll))
}

NR_root <- function(f, df, x, sq, w, i = 0, iter_max = 10000, tol = sqrt(2.22044604925e-16)){
    while(abs(f(x, sq, w)) > tol){
        x <- x - f(x, sq, w) / df(x, sq, w)
        i <- i + 1
        if (i == iter_max){break}
    }
    return(x)
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

sqrt_ginv <- function(X, tol = sqrt(2.22044604925e-16)){
    Xsvd <- svd(X);
    u <- Xsvd$u; s <- Xsvd$d; vh <- t(Xsvd$v)
    pos <- s > max(tol * s[1], 0)
    if (all(pos)){
        res <- t(vh) %*% diag(1 / s)**0.5 %*% t(u)
    } else if (all(!pos)){
        res <- matrix(0, nrow = nrow(X), ncol = ncol(X))
    } else {
        res <- t(vh)[, pos] %*% sqrt(diag(1 / s[pos])) %*% t(u[, pos])
    }
    return(res)
}
