# I expect that reference_quantiles is in sample_to_rescale as column
#' @export
scale_column_minmaxQ <- function(sample_to_rescale,
                                 scaling_values,
                                 colX,
                                 subtract_bg = TRUE) {
    scaling_values <- unlist(scaling_values, use.names = FALSE) # ensure unlisted.
    if (!is.list(scaling_values)) {
        # If the scaling_values are not a list, then it's the usual case of a vector
        # with one or two values, and classical relativisation is applied.
        sample_to_rescale <- scale_column_minmax(
            sample_to_rescale = sample_to_rescale,
            scaling_values = scaling_values,
            colX = colX,
            subtract_bg = subtract_bg
        )
        invisible(sample_to_rescale)
    }
    # If it's still a list, then it must be named, and each element might be a numeric or a function.
    scaling_values_funs <- lapply(scaling_values, function(x) {
        if (is.function(x)) {
            return(x)
        } else {
            return(function(y) {
                x
            })
        }
    })
    sample_to_rescale[, correction_low := ifelse(
        subtract_bg,
        scaling_values_funs[["negative"]](reference_quantiles),
        0
    )]
    sample_to_rescale[, correction_high := scaling_values_funs[["positive"]](reference_quantiles)]
    sample_to_rescale[, (colX) := (get(colX) - correction_low) / (correction_high - correction_low)]
    sample_to_rescale[, c("correction_low", "correction_high") := NULL]
    invisible(sample_to_rescale)
}
