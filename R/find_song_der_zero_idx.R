find_song_der_zero_idx <- function(all_info_lst, no_data, is.print = FALSE) {
    song_der_zero_idx_lst <- list()
    for(i in 1:no_data) {
        der_zero_idx <- which(diff(sign(all_info_lst[[i]][, 3]))!= 0)
        if (is.print) print(der_zero_idx)
        der_zero_idx_adj <- rep(0, length(der_zero_idx))

        ## choose beta1 closet to zero
        for (j in 1:length(der_zero_idx)) {

            if(abs(all_info_lst[[i]][der_zero_idx[j], 3]) > abs(all_info_lst[[i]][der_zero_idx[j]+1, 3])) {
                der_zero_idx_adj[j] <- der_zero_idx[j] + 1
            } else {
                der_zero_idx_adj[j] <- der_zero_idx[j]
            }
        }
        ## no end points
        end_idx <- der_zero_idx_adj %in% c(1, length(all_info_lst[[i]][, 3]))
        song_der_zero_idx_lst[[i]] <- der_zero_idx_adj[!end_idx]
    }
    song_der_zero_idx_lst
}
