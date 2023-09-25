## Song (2006) functions
Curve <- function(x, y, beta, CIp, imp) {
    # x <- data[ , 1]
    # y <- data[ , 2]
    l <- length(x)
    ma <- max(x)
    beta0 <- beta[ , 1]
    beta1 <- beta[ , 2]
    upperCI <- CIp[ , 1]
    lowerCI <- CIp[ , 2]
    # par(mfrow = c( 2, 1))
    ## nonparametric fitted curve with 95% CI
    plot(x, y, xlim = c(0, ma), type = 'p', ylim = c(-.9, 0.9),
         xlab = 'x', ylab = 'y', pch = 19)
    lines(x, beta0, lty = 4, lwd = 3, col = "green4")
    # title('Nonparametric Fitted Curve')
    pos <-as.vector(imp$position)
    vx <-as.vector(imp$x)
    ll <- length(vx)
    vx1 <-as.vector(imp$lCI)
    vx2 <-as.vector(imp$uCI)
    vy <- beta0[pos]

    for (i in 1:ll){
        x1 <- vx1[i]
        x2 <- vx2[i]
        y1 <- vy[i]
        polygon(c(x1, x1, x2, x2), c(0, 0.02, 0.02, 0) + 0.02*i, col = i,
                border = "white")
        points(vx[i], 0.02*i, pch = 18, col = i, cex = 1.5)
    }

    # ##first-derivative curve with 95% CI
    # plot(x, beta1, type = 'l', xlim = c(0, ma), xlab = 'x',
    #      ylab = 'First derivative')
    # title( 'Estimated First Derivative Curve with 95% Confidence Interval')
    # for (i in 1:( l - 1)){
    #     x1 <- x[i]
    #     x2 <- x[i + 1]
    #     y1 <- lowerCI[i]
    #     y2 <- lowerCI[i + 1]
    #     y3 <- upperCI[i]
    #     y4 <- upperCI[i + 1]
    #     polygon(c(x1, x1, x2, x2), c(y1, y3, y2, y4), col = "lightgrey")
    # }
    # cross <- rep(0, ll)
    # par(new = T)
    # plot(vx, cross, pch = "x", xlim = c( 0, ma), ylim = c( -2, 2),
    #      xlab = '', ylab = '')
    # abline(0)
    # par(new = T)
    # plot(x, beta1, type = 'l', xlim = c( 0, ma), ylim = c( -2, 2),
    #      xlab = 'Chromosomal coordinates (in kb)', ylab = 'First derivative'
}
