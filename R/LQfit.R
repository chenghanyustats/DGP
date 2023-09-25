## Song (2006) functions
LQfit <- function(x, y, h) {
    # x <- data[ , 1]
    # y <- data[ , 2]
    l <- length(x)
    beta0 <- rep(0, l)
    beta1 <- rep(0, l)
    for(i in 1:l) {
        x.reg <- x - x[i]
        w <- dnorm(x - x[i], 0, h)
        fit <- lm(y ~ x.reg + I(x.reg ^ 2), weights = w)
        beta0[i] <- fit$coe[1]
        beta1[i] <- fit$coe[2]
    }
    beta <- cbind(beta0, beta1)
    return(beta)
}
