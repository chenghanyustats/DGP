## smooth taut string method
smqreg <- function (y, thr.const = 2.5, verbose = FALSE, bandwidth = -1, 
                    sigma = -1, localsqueezing = TRUE, squeezing.factor = 0.5, 
                    DYADIC = TRUE, firstlambda = 100, smqeps = 1/length(y), 
                    fsign = double(0), 
                    gensign = TRUE, tolerance = 1e-12, ...) {
    n <- length(y)
    lambda <- rep(firstlambda, n - 1)
    sigma <- mad((y[-1] - y[-n])/sqrt(2))
    if (verbose) 
        print(c("sigma is ", sigma))
    if (gensign) 
        fsign <- gensign(y, thr.const = thr.const, extrema.mean = TRUE, 
                         sigma = sigma, localsqueezing = localsqueezing, squeezing.factor = squeezing.factor)
    repeat {
        f <- smqnew(y = y, lambda = lambda, eps = smqeps, fsign = fsign, 
                    tolerance = tolerance)
        if (bandwidth < 0) {
            residuals <- y - f
            residuals <- residuals - mean(residuals)
            if (DYADIC) 
                residuals.wr <- multiwdwr(residuals, sqrt(thr.const * 
                                                              log(n)) * sigma)
            else residuals.wr <- nondymwdr(residuals, sqrt(thr.const * 
                                                               log(n)) * sigma)
        }
        if (verbose) {
            par(mfrow = c(2, 1))
            plot(y, col = "grey", ...)
            lines(f, col = "red")
            lines(residuals.wr, type = "l", col = "green")
            lines(fsign - 2, col = "blue")
            plot(lambda, ty = "b", ...)
            print("Press Enter")
            dum <- readline()
        }
        if (bandwidth > 0) 
            break
        ind <- (abs(residuals.wr) > 1e-10)
        ind2 <- ind[-1] | ind[-n]
        if (sum(ind) == 0) 
            break
        if (localsqueezing) 
            lambda[ind2] <- lambda[ind2] * squeezing.factor
        else lambda <- lambda * squeezing.factor
        if (min(lambda) < 1e-08) {
            par(mfrow = c(2, 1))
            plot(y, col = "grey", ...)
            lines(f, col = "red")
            lines(residuals.wr, type = "l", col = "green")
            lines(fsign - 2, col = "blue")
            plot(lambda, ty = "b", ...)
            print("ERROR, This should not happen")
            break
        }
    }
    list(y = f, nmax = findmod(f)$mod, loc = findmod(f)$x, sigma = sigma)
}
