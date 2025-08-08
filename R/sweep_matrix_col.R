#' Apply a Function to Each Column Using Column-wise Statistics
#'
#' Applies a function (e.g., subtraction) between each column of a matrix `x` and a
#' corresponding statistic value provided in `stats`.
#' This is equivalent to applying `sweep()` with `MARGIN = 2`, but faster.
#'
#' @param x
#' A numeric matrix. Each column will be modified using the corresponding value in `stats`.
#' @param stats
#' A numeric vector with length equal to `ncol(x)`. These are the per-column statistics to apply with `fun`.
#' @param fun
#' A function or a string naming a function (e.g., `"-"`, `"+"`, `"*"`, `"/"`).
#' Defaults to subtraction (`"-"`).
#'
#' @return
#' A numeric matrix of the same dimensions as `x`, with `fun(x[, i], stats[i])` applied column-wise.
#'
#' @examples
#' mat <- matrix(1:9, nrow = 3)
#' stats <- c(1, 2, 3)
#' sweep_matrix_col(mat, stats, fun = "-")
#'
#' n <- 1e3
#' n_colums <- 1e2
#' n_rows <- n / n_colums
#' m <- matrix(rnorm(n_colums), nrow = n_rows, ncol = n_colums)
#' stats <- rnorm(n_colums)
#' # microbenchmark::microbenchmark(
#' #     v1 = sweep(m, 2, stats, "-"),
#' #     v2 = sweep_matrix_col(m, stats, "-"),
#' #     times = 100,
#' #     check = "equal"
#' # )
#' @export
sweep_matrix_col <- function(x, stats, fun = "-") {
    # Ensure the function argument is treated as a function object
    fun <- match.fun(fun)
    # Replicate the stats vector to form a matrix with same dimensions as x
    stats_matrix <- matrix(stats, nrow = nrow(x), ncol = length(stats), byrow = TRUE)
    # Apply the operation element-wise
    fun(x, stats_matrix)
}
