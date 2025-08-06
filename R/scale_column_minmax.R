#' Min-Max Scaling for a Single Column
#'
#' This function scales a specified column of a data.table using provided
#' scaling values, neither usually nor necessarily to `[0, 1]`.
#' If a single value is provided, then only one median given which should be centered around 0.
#' The transformation is performed in-place.
#'
#' @param sample_to_rescale A `data.table` containing the data to be rescaled.
#' @param scaling_values A numeric vector of length 1 or 2:
#' - If length 2: used as (min, max) for scaling.
#' - If length 1: interpreted as median and centered by subtracting it (see below).
#' @param colX The name of the column (as a string) to scale.
#' @param subtract_bg Logical. If TRUE, the minimum value (`low`) is subtracted from the data.
#' Defaults to TRUE. Ignored when centering on the median.
#'
#' @details
#' If `length(scaling_values) == 1`, the column is centered using:
#'   \deqn{x' = x - scaling\_values}
#' This is equivalent to subtracting the median.
#'
#' If `length(scaling_values) == 2`, then traditional min-max scaling is applied:
#'   \deqn{x' = (x - low) / (high - low)}
#'
#' The operation modifies the input `sample_to_rescale` in place.
#'
#' @return Invisibly returns the modified `sample_to_rescale`.
#'
#' @export
#' @keywords relativisation
scale_column_minmax <- function(sample_to_rescale,
                                scaling_values,
                                colX,
                                subtract_bg = TRUE) {
    # Careful, this function works always inplace!
    # Determine if we're centering instead of scaling
    mediancenter <- FALSE
    if (length(scaling_values) == 1) {
        # Then only one median given which should be centered around 0
        #
        #
        # For why the following works, see
        # ccc.transforms.transforms_functional.py scale_low_high_01():
        #       If there is only one value (len(x) == 1), then
        #         low -> x
        #         high -> low + 1
        #       which results in
        #         (tensor - low) / (high - low) =
        #         (tensor - low) / (low + 1 - low) =
        #         (tensor - low) / (1) = tensor - low
        scaling_values <- c(scaling_values, scaling_values + 1)
        mediancenter <- TRUE
    } else if (length(scaling_values) > 2) {
        stop("Expected 1 or 2 scaling values, got more than 2.")
    }

    low <- min(scaling_values)
    high <- max(scaling_values)

    # Subtract low value unless explicitly disabled
    low_subtractbg <- ifelse(subtract_bg || mediancenter, low, 0)

    # In-place scaling of the column
    data.table::set(
        sample_to_rescale,
        i = NULL,
        j = colX,
        value = (sample_to_rescale[, colX, with = FALSE] - low_subtractbg) / (high - low)
    )

    invisible(sample_to_rescale)
}
