#' Relative Scaling for a Single Column
#'
#' Scales a column of a `data.table` by subtracting a background value and dividing by a specified scale.
#' This method is commonly used for arcsinh-transformed cytometry data.
#' The operation is performed in-place.
#'
#' @param sample_to_rescale A `data.table` containing the data to be scaled.
#' @param scaling_values A numeric vector:
#' - If of length 1: interpreted as the background (e.g., median of negative population), and the column is centered by subtracting it and dividing by 1.
#' - If of length 2: interpreted as (low, high) values for subtraction and division.
#' @param colX The name (character) of the column in `sample_to_rescale` to scale.
#' @param subtract_bg Logical. If `TRUE` (default), the background (low value) is subtracted. If `FALSE`, it is assumed to be 0.
#'
#' @details
#' The formula used is:
#' \deqn{x' = (x - low) / high}
#' If `subtract_bg = FALSE`, then \eqn{low = 0}.
#'
#' For single-valued `scaling_values`, this behaves like centering the data around the median.
#'
#' @return Invisibly returns the modified `sample_to_rescale`.
#' @export
#' @keywords relativisation
scale_column_relative <- function(sample_to_rescale,
                                  scaling_values,
                                  colX,
                                  subtract_bg = TRUE) {
    # Careful, this function works always inplace!

    if (length(scaling_values) == 1) {
        # Then only one median given which should be centered around 0
        low <- scaling_values
        high <- 1
    } else if (length(scaling_values) == 2) {
        low <- min(scaling_values)
        high <- max(scaling_values)
    } else {
        stop("Expected 1 or 2 scaling values; got more than 2.")
    }

    # Optionally set background to zero
    if (!subtract_bg) {
        low <- 0
    }

    # Apply transformation
    data.table::set(
        sample_to_rescale,
        i = NULL,
        j = colX,
        value = (sample_to_rescale[, colX, with = FALSE] - low) / (high)
    )

    invisible(sample_to_rescale)
}
