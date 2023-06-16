#' @export 
rescale <- function(sample_to_rescale,
                    extracted_mfi,
                    feature_unified_dict,
                    missing_feature = c("minmax", "center_median"),
                    inplace_datatable = FALSE,
                    known_missing_features = c("FS INT", "FS TOF", "SS INT"),
                    scale_column_fun = scale_column_relative, ...) {
    if (!inplace_datatable || !data.table::is.data.table(sample_to_rescale)) {
        # That here effectively generates a new datatable
        # If the matrix _IS_ already a datatable, effectively a copy is created.
        # If the matrix is NOT created, this is _always_ called
        sample_to_rescale <- data.table::as.data.table(sample_to_rescale)
    }

    given_scale_column_fun <- scale_column_fun

    for (colX in colnames(sample_to_rescale)) {
        scale_column_fun <- given_scale_column_fun
        if (!colX %in% rownames(feature_unified_dict)) {
            if (missing_feature[1] == "minmax") {
                extracted_values <- c(
                    sample_to_rescale[, lapply(.SD, max), .SDcols = colX],
                    sample_to_rescale[, lapply(.SD, min), .SDcols = colX]
                )
                if (!colX %in% known_missing_features) {
                    warning(paste0(
                        "No MFI feature values for \'", colX, "\', used min and max of the data for standardization"
                    ))
                }
                scale_column_fun <- scale_column_minmax
            } else if (missing_feature[1] == "center_median") {
                extracted_values <- sample_to_rescale[, lapply(.SD, median), .SDcols = colX]
                if (!colX %in% known_missing_features) {
                    warning(paste0(
                        "No MFI feature values for \'", colX, "\', used median of the data for centering at that median."
                    ))
                }
                scale_column_fun <- scale_column_minmax
            } else {
                stop("Missing '", colX, "' rescaling")
            }
            extracted_values <- unlist(extracted_values)
        } else {
            extracted_values <- unlist(extracted_mfi[, feature_unified_dict[colX, "unified_single_staining"]])
        }

        scale_column_fun(
            sample_to_rescale = sample_to_rescale,
            scaling_values = extracted_values,
            colX = colX,
            ...
        )
    }
    return(sample_to_rescale)
}
