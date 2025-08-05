#' Rescale Sample Using Extracted MFIs
#'
#' Rescales all features (columns) in a cytometry sample using median fluorescence intensities (MFIs) extracted from controls.
#'
#' @inheritParams extract_marker_mfi_list
#' @inheritParams rescale_named
#'
#' @export
rescale_extracted <- function(sample_to_rescale,
                              extracted_mfi,
                              column_negative = c("negative", "unstained"),
                              column_positive = "positive",
                              missing_feature = c("minmax", "center_median"),
                              inplace_datatable = FALSE,
                              scale_column_fun = scale_column_minmax, ...) {
    # Convert extracted MFI table into named list of feature -> MFI (negative/positive)
    mfis_namedlist <- extract_marker_mfi_list(
        sample_to_rescale = sample_to_rescale,
        extracted_mfi = extracted_mfi,
        column_negative = column_negative,
        column_positive = column_positive
    )

    return(
        # Apply the actual rescaling using the named MFI list
        rescale_named(
            sample_to_rescale = sample_to_rescale,
            extracted_mfi_namedlist = mfis_namedlist,
            missing_feature = missing_feature,
            inplace_datatable = inplace_datatable,
            scale_column_fun = scale_column_fun,
            ...
        )
    )
}
