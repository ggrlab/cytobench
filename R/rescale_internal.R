rescale_internal <- function(sample_to_rescale,
                             extracted_mfi,
                             feature_unified_dict,
                             inplace_datatable = FALSE,
                             dt_fun) {
    if (!inplace_datatable || !data.table::is.data.table(sample_to_rescale)) {
        # That here effectively generates a new datatable
        # If the matrix _IS_ already a datatable, effectively a copy is created.
        # If the matrix is NOT created, this is _always_ called
        sample_to_rescale <- data.table::as.data.table(sample_to_rescale)
    }

    for (colX in colnames(sample_to_rescale)) {
        if (!colX %in% rownames(feature_unified_dict)) {
            stop("No MFI feature values for \'", colX, "\'")
        } else {
            extracted_values <- unlist(extracted_mfi[, feature_unified_dict[colX, "unified_single_staining"]])
        }

        dt_fun(sample_to_rescale, colX, extracted_values)
    }
    return(sample_to_rescale)
}
