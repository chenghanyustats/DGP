## log_post_t_thoery_matern

log_post_t_theory_matern <- function(t, y, x, Kff, A, lambda, l, sig2, nu = 1.5,
                                     shape1 = 1, shape2 = 1, a = 0, b = 2) {
    n <- length(x)
    tau2 <- sig2 / (n * lambda)
    Kdf <- computeCovDer1(idx1 = t, idx2 = x, ker_der1_fcn = matern_ker_der1,
                          tau = 1, l = l, nu = nu)
    if (nu == 1.5) {
        Kdd <- 3 / (l ^ 2)
    } else {
        Kdd <- 5 / (3 * l ^ 2)
    }

    a_mat <- t(Kdf) * (1 / sqrt(Kdd))
    mu_t <- t(a_mat) %*% solve(A, y) * sqrt(Kdd)
    dd <- 1 - quad.form.inv(A, a_mat)
    sig2_t <- tau2 * dd * Kdd
    cholA <- chol(A)
    log_det_A <- 2 * sum(log(diag(cholA)))

    # log_C <- (-n / 2) * log(2 * pi * tau2) -
    #     (1 / 2) * log(det(A)) - quad.form.inv(A, y) / (2 * tau2) +
    #     (1 / 2) * (log(tau2) + log(Kdd))

    log_C <- (-1 / 2) * (n * log(2 * pi * tau2) + log_det_A +
                             quad.form.inv(A, y) / tau2 -
                             (log(tau2) + log(Kdd)))

    log_den_t <- (-1 / 2) * (log(sig2_t) + crossprod(mu_t) / sig2_t)
    log_prior <- log_dgbeta(t, shape1, shape2, a, b)

    return(log_C + log_den_t + log_prior)
}
