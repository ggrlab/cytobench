#' Rescale Sample Using Extracted MFIs
#'
#' Rescales all features (columns) in a cytometry sample using median fluorescence intensities (MFIs) extracted from controls.
#'
#' @inheritParams extract_marker_mfi_list
#' @inheritParams rescale_named
#'
#' @export
#' @keywords relativisation
#' @examples
#' set.seed(42)
#' fs_ss <- simulate_cd3()
#' tmpdir <- local_tempdir_time()
#' flowCore::write.flowSet(fs_ss, tmpdir)
#'
#' extracted_mfis_singlestain <- extract_mfi(
#'     tmpdir,
#'     regex_singlestain = "(-(CD3-.*)|(none))\\.fcs$"
#' )
#'
#' rescaled_sample <- rescale_extracted(
#'     sample_to_rescale = flowCore::exprs(fs_ss[["sample0_12-panel"]]),
#'     extracted_mfi = extracted_mfis_singlestain
#' )
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
