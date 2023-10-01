se_prime <- function(x1, x2, h){
    n1 <- length(x1)
    n2 <- length(x2)
    XX <- matrix(rep(x2, n1), ncol = n1) - t(matrix(rep(x1, n2), ncol = n2))
    return(-2 * XX/h^2 * matrix(exp(-XX^2/h^2), ncol = n1))
}
