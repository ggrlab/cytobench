#' Extract Marker MFIs as Named List for Rescaling
#'
#' This function extracts median fluorescence intensities (MFIs) for specified features
#' from a tibble, returning a named list that maps each feature (e.g., "FITC-A") to a
#' named numeric vector of its negative and positive control values. If a feature
#' present in the sample is missing from the MFI table, it is added with `NA`.
#'
#' @param sample_to_rescale A `flowFrame` or similar object with marker columns to rescale.
#' @param extracted_mfi A tibble with columns: `feature`, `negative`, `positive` (or custom).
#' @param column_negative
#' Character vector of acceptable column names for the negative control.
#' Defaults to `c("negative", "unstained")`.
#' @param column_positive Name of the column for the positive control. Defaults to `"positive"`.
#'
#' @return A named list where each element corresponds to a marker (column) in `sample_to_rescale`.
#'         Each list element is a numeric vector of length 2 with names `negative` and `positive`.
#'         If a column in `sample_to_rescale` has no entry in `extracted_mfi`, it is filled with `NA`.
#'
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
            return(unlist(current_feature_mfi[, c(column_negative[1], column_positive[1])]))
        }
    )

    # Identify markers in the sample that are missing from the MFI table
    cn_sample <- flowCore::colnames(sample_to_rescale)
    cn_sample_missing <- cn_sample[!cn_sample %in% extracted_mfi[["feature"]]]

    # Add missing features with NA entries
    mfis_namedlist <- c(
        extracted_mfi_namedlist,
        setNames(rep(list(NA), length(cn_sample_missing)), cn_sample_missing)
    )

    return(mfis_namedlist)
}
