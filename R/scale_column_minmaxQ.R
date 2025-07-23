#' Rescale a Column Using Quantile-Based Min-Max Normalization
#'
#' This function rescales a column in a data.table using min-max scaling, where the
#' minimum and maximum are either fixed numeric values or functions that return values
#' based on `reference_quantiles`, which must be present in the `sample_to_rescale`.
#' It supports standard two-point scaling and more flexible quantile-based corrections.
#'
#' @param sample_to_rescale A data.table containing the data to be rescaled. Must include `reference_quantiles`.
#' @param scaling_values Either a numeric vector of length 2 (`c(min, max)`), or a named list
#'        with elements `"negative"` and `"positive"`, each being a numeric value or a function
#'        returning a value when applied to `reference_quantiles`.
#'        The scaling function should take a quantile value and return the corresponding high value.
#'
#' @param colX Character string specifying the column in `sample_to_rescale` to be rescaled.
#' @param subtract_bg Logical; if `TRUE`, the lower bound (negative control) is subtracted before scaling.
#'
#' @return Invisibly returns the modified `sample_to_rescale`, with the column `colX` scaled.
#' @export
scale_column_minmaxQ <- function(sample_to_rescale,
                                 scaling_values,
                                 colX,
                                 subtract_bg = TRUE) {
    # Flatten input in case it's a list with single numeric vector
    scaling_values <- unlist(scaling_values, use.names = FALSE)

    # Case 1: Standard min-max scaling (numeric vector input)
    if (!is.list(scaling_values)) {
        # Use default relativisation logic with static min and max values
        sample_to_rescale <- scale_column_minmax(
            sample_to_rescale = sample_to_rescale,
            scaling_values = scaling_values,
            colX = colX,
            subtract_bg = subtract_bg
        )
        return(invisible(sample_to_rescale))
    }

    # Case 2: scaling_values is a named list, e.g., list(negative = ..., positive = ...)
    # Each entry can be a function or a static value.
    # Convert all entries into functions: values become constant-returning functions.
    scaling_values_funs <- lapply(scaling_values, function(x) {
        if (is.function(x)) {
            return(x)
        } else {
            return(function(y) x) # Return a function that always returns x
        }
    })

    # Apply the lower and upper scaling functions to `reference_quantiles`
    # This assumes `reference_quantiles` is a column in `sample_to_rescale`
    # with the quantile values of the reference marker.
    sample_to_rescale[, correction_low := ifelse(
        subtract_bg,
        scaling_values_funs[["negative"]](reference_quantiles),
        0 # No background subtraction if FALSE
    )]
    # Using the reference quantiles, get the respective high value from the scaling function.
    # The scaling function should take a quantile value and return the corresponding high value.
    sample_to_rescale[, correction_high := scaling_values_funs[["positive"]](reference_quantiles)]

    # Rescale target column using dynamic low and high corrections
    sample_to_rescale[, (colX) := (get(colX) - correction_low) / (correction_high - correction_low)]

    # Clean up helper columns
    sample_to_rescale[, c("correction_low", "correction_high") := NULL]

    return(invisible(sample_to_rescale))
}
