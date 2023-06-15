rescale_internal <- function(sample_to_rescale,
                           #   Matrix, the actual sample
                           extracted_mfi,
                           # Result of single_staining_extract_mfi, e.g.
                           # extracted_mfi <- single_staining_extract_mfi(
                           # single_stainings = read_files_asinh[
                           #     sprintf("Align_01 CD3 Lot_A _EX Tube %02d %03d LC.csv", 1:10, 1:10)
                           # ],
                           # sample_unstained = read_files_asinh[[
                           #     sprintf("Align_01 CD3 Lot_A _EX Tube %02d %03d LC.csv", 11, 11)
                           # ]]
                           feature_unified_dict,
                           missing_feature = c("minmax", "center_median"),
                           inplace_datatable = FALSE,
                           known_missing_features = c("FS INT", "FS TOF", "SS INT")) {
    if (!inplace_datatable || !data.table::is.data.table(sample_to_rescale)) {
        # That here effectively generates a new datatable
        # If the matrix _IS_ already a datatable, effectively a copy is created.
        # If the matrix is NOT created, this is _always_ called
        sample_to_rescale <- data.table::as.data.table(sample_to_rescale)
    }

    for (colX in colnames(sample_to_rescale)) {
        if (!colX %in% rownames(feature_unified_dict)) {
            if (missing_feature[1] == "minmax") {
                both_values <- c(
                    sample_to_rescale[, lapply(.SD, max), .SDcols = colX],
                    sample_to_rescale[, lapply(.SD, min), .SDcols = colX]
                )
                if (!colX %in% known_missing_features) {
                    warning(paste0(
                        "No MFI feature values for \'", colX, "\', used min and max of the data for minmax-standardization"
                    ))
                }
            } else if (missing_feature[1] == "center_median") {
                both_values <- sample_to_rescale[, lapply(.SD, median), .SDcols = colX]
                if (!colX %in% known_missing_features) {
                    warning(paste0(
                        "No MFI feature values for \'", colX, "\', used median of the data for centering at that median."
                    ))
                }
            } else {
                stop("Missing '", colX, "' rescaling")
            }
            both_values <- unlist(both_values)
        } else {
            both_values <- unlist(extracted_mfi[, feature_unified_dict[colX, "unified_single_staining"]])
        }
        if (length(both_values) == 1) {
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
            both_values <- c(both_values, both_values + 1)
        } else if (length(both_values) > 2) {
            stop("Why are there more than 2 parameters?")
        }

        low <- min(both_values)
        high <- max(both_values)

        data.table::set(
            sample_to_rescale,
            i = NULL,
            j = colX,
            value = (sample_to_rescale[, colX, with = FALSE] - low) / (high - low)
        )
    }
    return(sample_to_rescale)
}
