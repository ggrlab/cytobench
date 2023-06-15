rescale_bg_relative <- function(sample_to_rescale,
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
                                known_missing_features = c("FS INT", "FS TOF", "SS INT"),
                                subtract_bg = TRUE) {
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
                        "No MFI feature values for \'", colX, "\', used min and max of the data for (x-min(single_staining))/max(single_staining)-standardization"
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
            low <- both_values
            high <- 1
        } else if (length(both_values) > 2) {
            stop("Why are there more than 2 parameters?")
        } else {
            low <- min(both_values)
            high <- max(both_values)
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
    }
    return(sample_to_rescale)
}
