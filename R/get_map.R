## get maximum a posteriori
get_map <- function(post_den, grid_t, hpdi) {
    # no_data <- nrow(post_den)
    grid_len <- length(grid_t)
    num_map <- hpdi$no_cluster
    map_est <- rep(0, num_map)
    for (k in 1:num_map) {
        lower <- which.min(abs(grid_t - hpdi$ci_lower[k]))
        upper <- which.min(abs(grid_t - hpdi$ci_upper[k]))
        map_est[k] <- grid_t[lower + which.max(post_den[lower:upper]) - 1]
    }
    map_est
}
