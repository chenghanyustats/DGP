
get_t <- function(b, a) {
    for (j in 1:length(b)) {
        if (j == 1) {
            idx <- which.min(abs(b[j] - a))
            b[idx] <- b[j] 
            idx_old <- idx
        } else {
            idx <- which.min(abs(b[j] - a))
            if (idx_old == idx) {
                if (abs(b[j] - a[idx_old]) < abs(b[idx_old] - a[idx_old])) {
                    b[idx] <- b[j] 
                    idx_old <- idx
                } else {
                    next
                }
            } else {
                b[idx] <- b[j] 
                idx_old <- idx
            }
        }
    }
    return(b[1:length(a)])
}
