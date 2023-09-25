## get hpd interval from density 
get_hpd_interval_from_den <- function(post_prob, grid_t,
                                      max_iter = 1000, tol = 0.01,
                                      step_size = 0.001,
                                      target_prob = 0.95) {
    sum_density <- sum(post_prob)
    max_density <- max(post_prob)
    current_density <- max_density
    prob_value <- 0
    count <- 1
    # sample_size <- length(samples)
    while(abs(prob_value - target_prob) > tol) {
        if (prob_value > target_prob) {
            current_density <- current_density + step_size
        } else {
            current_density <- current_density / sqrt(1 + count)
        }
        
        high_density_idx <- which(post_prob > current_density)
        clu_idx <- which(diff(high_density_idx) > 2)
        (prob_value <- sum(post_prob[high_density_idx]) / sum_density)
        
        ci_lower_b <- c(grid_t[high_density_idx][1], 
                        grid_t[high_density_idx][clu_idx + 1])
        ci_upper_b <- c(grid_t[high_density_idx][clu_idx],
                        grid_t[high_density_idx][length(high_density_idx)])
        no_clu <- length(ci_lower_b)
        # samples[samples > ci_lower & samples < ci_upper]
        count <- count + 1
        # print(prob_value)
        if (count == max_iter) break
    }
    return(list(no_cluster = no_clu, 
                ci_lower = ci_lower_b, ci_upper = ci_upper_b, 
                prob_value = prob_value,
                den_value = current_density))
}
