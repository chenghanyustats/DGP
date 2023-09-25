# get_sig2_t_new

get_sig2_t_new <- function(sig2, lambda, tau2 = NULL, h, Kdf, Ainv) {
    if(!is.null(tau2)) {
        return(tau2 * (1 / h ^ 2 - quad.tform(Ainv, Kdf)))
    } else {
        sig2 / (nrow(A) * lambda) * (1 / h ^ 2 - quad.tform(Ainv, Kdf))
    }
}
