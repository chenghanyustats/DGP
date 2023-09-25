# get_mu_t_new

get_mu_t_new <- function(y, Kdf, Ainv) {
    Kdf %*% Ainv %*% y
}
