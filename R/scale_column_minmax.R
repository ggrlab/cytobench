#' @export 
scale_column_minmax <- function(sample_to_rescale,
                                scaling_values,
                                colX,
                                subtract_bg = TRUE) {
    # Careful, this function works always inplace!


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
    } else if (length(scaling_values) > 2) {
        stop("Why are there more than 2 parameters?")
    }

    low <- min(scaling_values)
    high <- max(scaling_values)

    data.table::set(
        sample_to_rescale,
        i = NULL,
        j = colX,
        value = (sample_to_rescale[, colX, with = FALSE] - low) / (high - low)
    )
    invisible(sample_to_rescale)
}
