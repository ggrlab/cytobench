#' Extract Marker MFIs as Named List for Rescaling
#'
#' This function extracts median fluorescence intensities (MFIs) for specified features
#' from a tibble, returning a named list that maps each feature (e.g., "FITC-A") to a
#' named numeric vector of its negative and positive control values. Features present
#' in `sample_to_rescale` but missing from `extracted_mfi` are included with `NA`.
#'
#' @param sample_to_rescale A `flowFrame` or other object with marker (column) names to be rescaled.
#' @param extracted_mfi A `tibble` (or coercible to one) with at least three columns:
#'   \describe{
#'     \item{`feature`}{Character. Names of the features (e.g., `"FITC-A"`).}
#'     \item{`negative`}{Numeric. MFI values for negative/unlabeled control.}
#'     \item{`positive`}{Numeric. MFI values for positive/labeled control.}
#'   }
#' @param column_negative Character vector. Column name as  negative control.
#'   Defaults to `c("negative", "unstained")`. First element is used.
#' @param column_positive Character. Column name used as positive control. Defaults to `"positive"`. First element is used.
#'
#' @return A named list. Each entry corresponds to a marker (column) in `sample_to_rescale`,
#'   containing a numeric vector of length two with names `"negative"` and `"positive"`.
#'   Missing markers are filled with `NA`.
extract_marker_mfi_list <- function(
    sample_to_rescale,
    extracted_mfi,
    column_negative = c("negative", "unstained"),
    column_positive = "positive") {
    feature <- NULL # to avoid R CMD check note about undefined global variable

    # Ensure `extracted_mfi` is a tibble to prevent data.table-specific issues
    extracted_mfi <- tibble::as_tibble(extracted_mfi)

    # Build named list of MFIs for each marker/feature
    extracted_mfi_namedlist <- sapply(
        as.character(extracted_mfi[["feature"]]),
        simplify = FALSE,
        function(feature_x) {
            current_feature_mfi <- extracted_mfi |>
                dplyr::filter(feature == feature_x)
            # Extract the two relevant MFI values and return as named vector
            unlist(current_feature_mfi[, c(column_negative[1], column_positive[1])])
        }
    )

    # Identify any marker in sample that is missing from MFI table
    cn_sample <- flowCore::colnames(sample_to_rescale)
    cn_sample_missing <- cn_sample[!cn_sample %in% extracted_mfi[["feature"]]]

    # Fill missing markers with NA
    mfis_namedlist <- c(
        extracted_mfi_namedlist,
        stats::setNames(rep(list(NA), length(cn_sample_missing)), cn_sample_missing)
    )

    return(mfis_namedlist)
}
