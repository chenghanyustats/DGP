nks_ci <- function(no_data, all_info_lst, song_der_zero_pt_lst) {
    ci_pts_song_lst <- list()
    n <- nrow(all_info_lst[[1]])
    for (i in 1:no_data) {
        # print(paste("i=", i))
        stationary_pt <- all_info_lst[[i]][song_der_zero_pt_lst[[i]], ]
        if(length(song_der_zero_pt_lst[[i]]) == 1) {
            stationary_pt_upper <- stationary_pt[4]
            stationary_pt_lower <- stationary_pt[5]
        } else {
            stationary_pt_upper <- stationary_pt[, 4]
            stationary_pt_lower <- stationary_pt[, 5]
        }
        
        upper_vec <- c()
        lower_vec <- c()
        for(j in 1:length(song_der_zero_pt_lst[[i]])) {
            idx <- song_der_zero_pt_lst[[i]][j]
            if (all_info_lst[[i]][idx, 3] - 
                all_info_lst[[i]][idx-1, 3] < 0) {
                
                
                upper_idx_search <- which(diff(rev(all_info_lst[[i]][1:(idx-1), 3])) < 0)
                
                if (length(upper_idx_search) != 0) {
                    upper_idx <- all_info_lst[[i]][(idx-1):(idx-1-upper_idx_search[1]), 3][all_info_lst[[i]][(idx-1):(idx-1-upper_idx_search[1]), 3] > 
                                                                                               stationary_pt_upper[j]][1]
                } else {
                    upper_idx <- rev(all_info_lst[[i]][1:(idx-1), 3][all_info_lst[[i]][1:(idx-1), 3] > 
                                                                         stationary_pt_upper[j]])[1]
                }
                
                if (!is.na(upper_idx)) {
                    lower_pt <- all_info_lst[[i]][which(all_info_lst[[i]][, 3] == upper_idx), 1]
                }
                
                
                
                lower_idx_search <-  which(diff(all_info_lst[[i]][(idx+1):n, 3]) > 0)
                
                
                if(length(lower_idx_search) != 0) {
                    lower_idx <- all_info_lst[[i]][(idx+1):(idx+lower_idx_search[1]), 3][all_info_lst[[i]][(idx+1):(idx+lower_idx_search[1]), 3] < 
                                                                                             stationary_pt_lower[j]][1]
                } else{
                    lower_idx <- all_info_lst[[i]][(idx+1):n, 3][all_info_lst[[i]][(idx+1):n, 3] < 
                                                                     stationary_pt_lower[j]][1]
                }
                
                
                if (!is.na(lower_idx)) {
                    upper_pt <- all_info_lst[[i]][which(all_info_lst[[i]][, 3] == lower_idx), 1]
                }
                
                
                if(is.na(upper_idx) & !is.na(lower_idx)) {
                    
                    lower_pt <- 2 * stationary_pt[j, 1] - upper_pt
                }
                
                if(!is.na(upper_idx) & is.na(lower_idx)) {
                    upper_pt <- 2 * stationary_pt[j, 1] - lower_pt
                }
                
                if(is.na(upper_idx) & is.na(lower_idx)) {
                    lower_pt <- NA
                    upper_pt <- NA
                }
                
                
            } else {
                upper_idx_search <- which(diff(all_info_lst[[i]][(idx+1):n, 3]) < 0)
                if (length(upper_idx_search) != 0) {
                    upper_idx <- all_info_lst[[i]][(idx+1):(idx+upper_idx_search[1]), 3][all_info_lst[[i]][(idx+1):(idx+upper_idx_search[1]), 3] > 
                                                                                             stationary_pt_upper[j]][1]
                } else {
                    upper_idx <- all_info_lst[[i]][(idx+1):n, 3][all_info_lst[[i]][(idx+1):n, 3] >
                                                                     stationary_pt_upper[j]][1]
                }
                
                if (!is.na(upper_idx)) {
                    upper_pt <- all_info_lst[[i]][which(all_info_lst[[i]][, 3] == upper_idx), 1]
                }
                
                lower_idx_search <- which(diff(rev(all_info_lst[[i]][1:(idx-1), 3])) > 0)
                if(length(lower_idx_search) != 0) {
                    lower_idx <- all_info_lst[[i]][(idx-1):(idx-lower_idx_search[1]), 3][all_info_lst[[i]][(idx-1):(idx-lower_idx_search[1]), 3] < 
                                                                                             stationary_pt_lower[j]][1]
                } else{
                    lower_idx <- rev(all_info_lst[[i]][1:(idx-1), 3][all_info_lst[[i]][1:(idx-1), 3] <
                                                                         stationary_pt_lower[j]])[1]
                }
                
                if (!is.na(lower_idx)) {
                    lower_pt <- all_info_lst[[i]][which(all_info_lst[[i]][, 3] == lower_idx), 1]
                }
                
                
                if(is.na(upper_idx) & !is.na(lower_idx)) {
                    upper_pt <- 2 * stationary_pt[j, 1] - lower_pt
                }
                
                if(!is.na(upper_idx) & is.na(lower_idx)) {
                    lower_pt <- 2 * stationary_pt[j, 1] - upper_pt
                }
                
                if(is.na(upper_idx) & is.na(lower_idx)) {
                    lower_pt <- NA
                    upper_pt <- NA
                }
            }
            upper_vec <- c(upper_vec, upper_pt)
            lower_vec <- c(lower_vec, lower_pt)
        }
        
        
        ci_pts_song_lst[[i]] <- cbind(lower_vec, upper_vec)
    }
    return(ci_pts_song_lst)
}
