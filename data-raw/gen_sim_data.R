################################################################################
# GP with Derivative Information Simulation Study: Data Generating             #                                                              #
# "./data/GenSimData.R"                                                        #
# Cheng-Han Yu                                                                 #
################################################################################
# library(here)
# source("data-raw/seed.R")

## data list
YY <- list()

## size
n <- 100  # 500 and 1000

## number of replicate
no_data <- 100

## standard deviation of error noise
sig <- .1


## input range
x_a <- 0
x_b <- 1


## true reg function
regfcn <- function(x) {
    # reg <- 0.5 * 6 ^ x * sin(4.1 / x)
    reg <-sqrt(x*(1-x))*sin(2 * pi / (x + 0.5))
    return(reg)
}


## first derivative
fn <- function(x) {
    # 0.5 * log(6) * sin(4.1 / x) * 6 ^ x - (4.1 / 2 * cos(4.1 / x) * 6 ^ (x) / (x ^ 2))
    (1-2*x)*sin(2*pi/(x+1/2)) / (2 * sqrt(x*(1-x))) - 2*pi*sqrt(x*(1-x))*cos(2*pi/(x+1/2)) / ((x+1/2)^2)
}

## stationary point
cri_pts <- c(uniroot(fn, c(x_a, 0.2))$root, uniroot(fn, c(0.2, 0.4))$root, uniroot(fn, c(0.4, x_b))$root)



## produce data
for (k in 1:no_data) {
    x <- seq(x_a, x_b, length.out = n)
    set.seed(seeds[k])
    y <- regfcn(x) + rnorm(n, 0, sig)
    YY[[k]] <- list(x = x, y = y)
}



## plotting
par(mfrow = c(1, 1), mar = c(4, 4, 2, 0))

idx <- seq(x_a, x_b, length.out = 500)
truefcn <- regfcn(idx)
plot(idx, truefcn, col = "red", type = "l", xlim = c(x_a - 0.05, x_b + 0.05),
     ylim = c(-.9, 0.9), lwd = 2, main = "Simulated Data of Size 100",
     las = 1, xlab = "x", ylab = "f(x)")
points(x, y, pch = 3)
points(cri_pts, regfcn(cri_pts), pch = "_", lwd = .1,
       col = "green", cex = 3)


# save(YY, sig, x_a, x_b, cri_pts, no_data,
#      file = "./data/sim_data_n100.RData")
