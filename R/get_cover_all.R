get_cover_all <- function(x, y, n, x_new, cri_pts) {
    powell <- Rsolnp::solnp(pars = c(0.5, 0.5), fun = log_post_lambda,
                            LB = c(0, 0), UB = c(1, 1),
                            x = x, y = y, n = n)
    h <- powell$pars[1]
    lambda <- powell$pars[2]
    return(get_cover(h, lambda, x, x_new, n, y, cri_pts))
}
