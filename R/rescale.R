#' Rescale Sample Using Unified Feature Dictionary
#'
#' This function rescales a sample (data.table) using extracted MFI values and a unified feature mapping.
#' If a feature is not listed in the dictionary, fallback methods such as min-max or median-centering are used.
#'
#' @param feature_unified_dict A `data.frame` where rownames are marker names in `sample_to_rescale`, and
#' one column (`unified_single_staining`) maps to column names in `extracted_mfi`.
#'
#' @inheritParams extract_marker_mfi_list
#' @inheritParams rescale_named
#' @param known_missing_features
#' A character vector of known features that may not have MFI values in `extracted_mfi`.
#' These features will not raise warnings when missing, but will use fallback methods.
#' Defaults to `c("FS INT", "FS TOF", "SS INT").
#'
#' @return A `data.table` with rescaled columns.
#' @export
rescale <- function(sample_to_rescale,
                    extracted_mfi,
                    feature_unified_dict,
                    missing_feature = c("minmax", "center_median"),
                    inplace_datatable = FALSE,
                    known_missing_features = c("FS INT", "FS TOF", "SS INT"),
                    scale_column_fun = scale_column_minmax,
                    ...) {
    # Ensure we work on a data.table (create a copy if needed)
    if (!inplace_datatable || !data.table::is.data.table(sample_to_rescale)) {
        sample_to_rescale <- data.table::as.data.table(sample_to_rescale)
    }

    # Preserve original function so it can be reassigned per-column if necessary
    given_scale_column_fun <- scale_column_fun

    for (colX in colnames(sample_to_rescale)) {
        scale_column_fun <- given_scale_column_fun

        # Determine whether this feature exists in the unified dictionary
        if (!colX %in% rownames(feature_unified_dict)) {
            # Feature not in dictionary: apply fallback method
            if (missing_feature[1] == "minmax") {
                extracted_values <- c(
                    sample_to_rescale[, lapply(.SD, max), .SDcols = colX],
                    sample_to_rescale[, lapply(.SD, min), .SDcols = colX]
                )
                if (!colX %in% known_missing_features) {
                    warning(paste0(
                        "No MFI feature values for '", colX,
                        "', used min and max of the data for standardization."
                    ))
                }
                scale_column_fun <- scale_column_minmax
            } else if (missing_feature[1] == "center_median") {
                extracted_values <- sample_to_rescale[, lapply(.SD, median), .SDcols = colX]
                if (!colX %in% known_missing_features) {
                    warning(paste0(
                        "No MFI feature values for '", colX,
                        "', used median of the data for centering."
                    ))
                }
                scale_column_fun <- scale_column_minmax
            } else {
                stop("No rescaling strategy defined for missing column: ", colX)
            }
        } else {
            # Extract rescaling values from MFI table using unified feature name
            unified_feature <- feature_unified_dict[colX, "unified_single_staining"]
            extracted_values <- extracted_mfi[, unified_feature]
        }

        extracted_values <- unlist(extracted_values)

        # Apply the rescaling function
        scale_column_fun(
            sample_to_rescale = sample_to_rescale,
            scaling_values = extracted_values,
            colX = colX,
            ...
        )
    }

    return(sample_to_rescale)
}
