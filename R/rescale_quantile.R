rescale_quantile <- function(
    sample_to_rescale,
    extracted_mfi,
    extracted_qfun,
    reference_marker = "AA750-A",
    column_negative = c("negative", "unstained"),
    column_positive = "positive",
    missing_feature = c("minmax", "center_median"),
    inplace_datatable = FALSE,
    scale_column_fun = scale_column_minmaxQ, ...) {
    mfis_namedlist <- extract_marker_mfi_list(
        sample_to_rescale = sample_to_rescale,
        extracted_mfi = extracted_mfi,
        column_negative = column_negative,
        column_positive = column_positive
    )
    mfis_namedlist <- lapply(mfis_namedlist, as.list)
    
    # For each feature in the extracted_qfun, we REPLACE the positive control value
    # with the quantile function value
    for(name_x in names(extracted_qfun)){
        mfis_namedlist[[name_x]][["positive"]] <- extracted_qfun[[name_x]]
    }
    return(
        rescale_named_quantiles(
            sample_to_rescale = sample_to_rescale,
            extracted_mfi_namedlist = mfis_namedlist,
            reference_marker = reference_marker,
            missing_feature = missing_feature,
            inplace_datatable = inplace_datatable,
            scale_column_fun = scale_column_fun,
            ...
        )
    )
}
