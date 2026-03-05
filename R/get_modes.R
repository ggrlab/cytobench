get_modes <- function(xy, force_1d = FALSE, ...) {
    if (ncol(xy) == 1) {
        return(hdrcde::hdr(xy, ...)$mode)
    } else if (ncol(xy) > 2) {
        stop("Only 1 or 2 dimensions are supported for mode calculation.")
    } else {
        if (force_1d) {
            modes <- apply(xy, 2, function(col_x) {
                res <- do.call(
                    hdrcde::hdr,
                    c(list(x = col_x), ...)
                )
                res$mode
            })
        } else {
            modes <- do.call(
                hdrcde::hdr.2d,
                c(
                    list(x = xy[, 1], y = xy[, 2]),
                    ...
                )
            )$mode
            names(modes) <- colnames(xy)
        }
    }
    return(modes)
}
