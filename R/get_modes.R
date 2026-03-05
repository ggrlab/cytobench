#' @title Compute Distribution Modes for 1D or 2D Data
#'
#' @description
#' Computes mode locations using highest-density-region estimators from
#' \\pkg{hdrcde} for one- or two-dimensional numeric data.
#'
#' @details
#' For one-dimensional input, this function calls \\code{hdrcde::hdr()}. For
#' two-dimensional input, it calls \\code{hdrcde::hdr.2d()} unless
#' \\code{force_1d = TRUE}, in which case each column is processed separately
#' with \\code{hdrcde::hdr()}. This is useful when marginal modes are preferred
#' over joint 2D mode estimates.
#'
#' @param xy A numeric vector, matrix, or data.frame with 1 or 2 columns.
#' @param force_1d Logical scalar; for 2D input, estimate one mode per column.
#' @param ... Additional arguments passed to \\code{hdrcde::hdr()} or
#' \\code{hdrcde::hdr.2d()}.
#'
#' @return
#' For 1D input, a numeric vector of mode location(s). For 2D input with
#' \\code{force_1d = FALSE}, a numeric vector of joint mode coordinates. For 2D
#' input with \\code{force_1d = TRUE}, either a named numeric vector (if one mode
#' per column) or a named list of numeric vectors (if a column has multiple
#' modes).
#'
#' @examples
#' if (requireNamespace("hdrcde", quietly = TRUE)) {
#'     set.seed(1)
#'     x <- c(rnorm(200, -2), rnorm(200, 2))
#'     get_modes(x)
#'     # Expected: mode estimates close to the two main peaks.
#'
#'     xy <- cbind(x = rnorm(400), y = rnorm(400, mean = 1))
#'     get_modes(xy)
#'     # Expected: one joint 2D mode coordinate pair.
#' }
#'
#' @keywords internal
get_modes <- function(xy, force_1d = FALSE, ...) {
    xy <- as.matrix(xy)
    if (!is.numeric(xy)) {
        stop("`xy` must contain numeric values.", call. = FALSE)
    }

    n_dims <- ncol(xy)
    if (!(n_dims %in% c(1L, 2L))) {
        stop("Only 1 or 2 dimensions are supported for mode calculation.",
            call. = FALSE
        )
    }

    if (n_dims == 1L) {
        return(do.call(hdrcde::hdr, c(list(x = xy[, 1]), list(...)))$mode)
    }
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
    return(modes)
}
