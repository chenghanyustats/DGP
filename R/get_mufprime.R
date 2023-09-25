## get muf
get_mufprime <- function(y, x, t, H0, tau, h, lambda) {
    n <- length(y)
    Kff_gp <- se_ker(H0 = H0, tau = tau, h = h)
    A_gp <- Kff_gp + diag((n * lambda), n)
    A_inv <- solve(A_gp)
    Kdf <- computeCovDer1(idx1 = t, idx2 = x, 
                          tau = tau, h = h) 
    A_inv_y <- A_inv %*% y
    mufprime <- Kdf %*% A_inv_y
}

