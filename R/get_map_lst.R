## get maximum a posteriori
get_map_lst <- function(post_den, grid_t, hpdi_1_2_lst) {
    no_data <- nrow(post_den)
    grid_len <- length(grid_t)

    map_est <- vector("list", length = no_data)
    for (i in 1:no_data) {
        num_map <- nrow(hpdi_1_2_lst[[i]])
        for (k in 1:num_map) {
            lower <- which.min(abs(grid_t - hpdi_1_2_lst[[i]][k, 1]))
            upper <- which.min(abs(grid_t - hpdi_1_2_lst[[i]][k, 2]))

            map_est[[i]][[k]] <- grid_t[lower + which.max(post_den[i, lower:upper]) - 1]
            # print(map_est[[i]][[k]])
        }

        # map_est[i, 1] <- grid_t[which.max(post_den[i, 1:313])]
        # map_est[i, 2] <- grid_t[which.max(post_den[i, 314:1000]) + 313]
    }
    map_est
}
