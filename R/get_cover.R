get_cover <- function(h, lambda, x, x_new, n, y, cri_pts){

    K_XX <- se(x, x, h)
    K_Xx <- se(x, x_new, h)
    K10_Xx <- se_prime(x, x_new, h)
    K11_xx <- se_2prime(x_new, x_new, h)

    eigendecomposition <- eigen(K_XX)
    eigenvector <- eigendecomposition$vector
    eigenvalue <- eigendecomposition$values

    # t0 <- c(0.0863302, 0.3095537, 0.7490745)
    t0 <- cri_pts

    f_hat <- as.vector(K_Xx %*% eigenvector %*% diag(1/(eigenvalue + n * lambda)) %*% t(eigenvector) %*% y)
    f_prime_hat <- as.vector(K10_Xx %*% eigenvector %*% diag(1/(eigenvalue + n * lambda)) %*% t(eigenvector) %*% y)
    sigma_hat <- sqrt(lambda * as.vector(y %*% eigenvector %*% diag(1/(eigenvalue + n * lambda)) %*% t(eigenvector) %*% y))
    cov_11 <- sigma_hat^2 / (n * lambda) * (K11_xx - K10_Xx %*% eigenvector %*% diag(1/(eigenvalue + n * lambda)) %*% t(eigenvector) %*% t(K10_Xx))
    cov_11 <- diag(cov_11)

    l <- sqrt(diag(K11_xx/cov_11)) * exp(-f_prime_hat^2/(2*cov_11))
    hpd <- get_hpd_interval_from_den(l, x_new)
    hpd_n <- hpd$no_cluster
    hpd_ci <- cbind(hpd$ci_lower, hpd$ci_upper)

    t_hat <- c()

    for (i in 1:hpd_n){
        x_can <- x_new[which(x_new == hpd_ci[i, 1]):which(x_new == hpd_ci[i, 2])]
        f_can <- f_hat[which(x_new == hpd_ci[i, 1]):which(x_new == hpd_ci[i, 2])]
        t_new <- sort(x_can[c(which(diff(sign(diff(f_can)))==-2)+1, which(diff(sign(diff(f_can)))==2)+1)])
        if (length(t_new) > 0){
            t_hat = c(t_hat, mean(t_new))
        }

    }

    K10_Xt <- se_prime(x, t_hat, h)
    K20_Xt <- se_2prime(x, t_hat, h)

    mu_hat <- as.vector(K_XX %*% eigenvector %*% diag(1/(eigenvalue + n * lambda)) %*% t(eigenvector) %*% y)
    delta <- K10_Xt %*% eigenvector %*% diag(1/(eigenvalue + n * lambda)) %*% t(eigenvector) %*% mu_hat
    f_prime2_hat <- as.vector(K20_Xt %*% eigenvector %*% diag(1/(eigenvalue + n * lambda)) %*% t(eigenvector) %*% mu_hat)
    bias <- delta / f_prime2_hat

    r <- sigma_hat * sqrt(diag(K10_Xt %*% eigenvector %*% diag(1/(eigenvalue + n * lambda)^2) %*% t(eigenvector) %*% t(K10_Xt))) / abs(f_prime2_hat)
    CI <- cbind(t_hat + bias + qnorm(0.05) * r, t_hat + bias + qnorm(0.95) * r)

    cov_01 <- rep(0, 3)
    for (i in 1:3){
        cov_01[i] <- t0[i] > CI[i,1] && t0[i] < CI[i,2]
    }

    CI <- cbind(t_hat + bias + qnorm(0.025) * r, t_hat + bias + qnorm(0.975) * r)
    cov_05 <- rep(0, 3)
    for (i in 1:3){
        cov_05[i] <- t0[i] > CI[i,1] && t0[i] < CI[i,2]
    }

    CI <- cbind(t_hat + bias + qnorm(0.005) * r, t_hat + bias + qnorm(0.995) * r)
    cov_001 <- rep(0, 3)
    for (i in 1:3){
        cov_001[i] <- t0[i] > CI[i,1] && t0[i] < CI[i,2]
    }

    bonferroni <- rep(0,3)
    CI <- cbind(t_hat + bias + qnorm(0.1/(2*length(t_hat))) * r, t_hat + bias + qnorm(1-0.1/(2*length(t_hat))) * r)
    bonferroni[1] <- any(cbind( CI[,1] <= t0[1] & CI[,2] >= t0[1] )) * any(cbind( CI[,1] <= t0[2] & CI[,2] >= t0[2] )) * any(cbind( CI[,1] <= t0[3] & CI[,2] >= t0[3] ))

    CI <- cbind(t_hat + bias + qnorm(0.05/(2*length(t_hat))) * r, t_hat + bias + qnorm(1-0.05/(2*length(t_hat))) * r)
    bonferroni[2] <- any(cbind( CI[,1] <= t0[1] & CI[,2] >= t0[1] )) * any(cbind( CI[,1] <= t0[2] & CI[,2] >= t0[2] )) * any(cbind( CI[,1] <= t0[3] & CI[,2] >= t0[3] ))

    CI <- cbind(t_hat + bias + qnorm(0.01/(2*length(t_hat))) * r, t_hat + bias + qnorm(1-0.01/(2*length(t_hat))) * r)
    bonferroni[3] <- any(cbind( CI[,1] <= t0[1] & CI[,2] >= t0[1] )) * any(cbind( CI[,1] <= t0[2] & CI[,2] >= t0[2] )) * any(cbind( CI[,1] <= t0[3] & CI[,2] >= t0[3] ))

    return(c("hpd_n" = hpd_n,
             "bonferroni_" = bonferroni,
             "cov_01_" = cov_01,
             "cov_05_" = cov_05,
             "cov_001_" = cov_001))
}
