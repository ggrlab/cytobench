#' @export 
scale_column_relative <- function(sample_to_rescale,
                                  scaling_values,
                                  colX,
                                  subtract_bg = TRUE) {
    # Careful, this function works always inplace!


    if (length(scaling_values) == 1) {
        # Then only one median given which should be centered around 0
        low <- scaling_values
        high <- 1
    } else if (length(scaling_values) > 2) {
        stop("Why are there more than 2 parameters?")
    } else {
        low <- min(scaling_values)
        high <- max(scaling_values)
    }

    if (!subtract_bg) {
        low <- 0
    }


    data.table::set(
        sample_to_rescale,
        i = NULL,
        j = colX,
        value = (sample_to_rescale[, colX, with = FALSE] - low) / (high)
    )
    invisible(sample_to_rescale)
}
