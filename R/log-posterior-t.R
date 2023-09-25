################################################################################
# Log posterior of t in the theory paper (given hyperparameters)
# with Beta prior
################################################################################
log_post_t <- function(t, y, x, Kff, A, sig2, lambda, h,
                            shape1 = 1, shape2 = 1, a = 0, b = 2) {

    Kdf <- computeCovDer1(t, x, ker_der1_fcn = se_ker_der1_no_tau, h = h)
    # B <- get_B(Kdf = Kdf, h = h)
    a_mat <- get_a_mat(Kdf = Kdf, h = h)
    at_Ainv <- t(a_mat) %*% chol2inv(A)
    mu_t <- get_mu_t(y, A, a_mat, h, at_Ainv = at_Ainv)
    sig2_t <-  get_sig2_t(sig2, lambda, A, a_mat, h, at_Ainv = at_Ainv)

    (-1 / 2) * (log(sig2_t) + mu_t ^ 2 / sig2_t) +
        log_dgbeta(t, shape1 = shape1, shape2 = shape2, a = a, b = b)
}



log_post_t_new <- function(t, y, x, Kff, A, Ainv, sig2, lambda, h, tau2,
                       shape1 = 1, shape2 = 1, a = 0, b = 2) {
    n <- length(y)
    Kdf <- computeCovDer1(t, x, ker_der1_fcn = se_ker_der1_no_tau, h = h)
    # B <- get_B(Kdf = Kdf, h = h)
    # a_mat <- get_a_mat(Kdf = Kdf, h = h)
    # at_Ainv <- t(a_mat) %*% chol2inv(A)
    mu_t <- get_mu_t_new(y, Kdf = Kdf, Ainv = Ainv)
    sig2_t <-  get_sig2_t_new(sig2, lambda, tau2, h, Kdf = Kdf, Ainv = Ainv)

    log_t <- (-1 / 2) * (log(sig2_t) + mu_t ^ 2 / sig2_t) +
        log_dgbeta(t, shape1 = shape1, shape2 = shape2, a = a, b = b)

    log_C <- (-n / 2) * log(2 * pi * tau2) -
        (1 / 2) * determinant(A, logarithm = TRUE)$modulus -
        quad.form(Ainv, y) - (2 * tau2) + (1 / 2) * log(tau2) - log(h)

    return(log_t + log_C)
}





