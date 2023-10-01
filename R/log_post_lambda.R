# Log marginal posterior likelihood
log_post_lambda = function(theta, x, y, n){
    K_XX <- se(x, x, theta[1])
    eigendecomposition <- eigen(K_XX)
    eigenvector <- eigendecomposition$vector
    eigenvalue <- eigendecomposition$values
    log_p <- n * log(1/n * as.vector(t(y) %*% eigenvector %*% diag(1/(eigenvalue + n * theta[2])) %*% t(eigenvector) %*% y)) + sum(log(eigenvalue + n * theta[2]))
    return(log_p)
}
