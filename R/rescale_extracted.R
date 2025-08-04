#' Rescale a sample using a named list of extracted MFIs
#'
#' This function rescales all columns of a sample using the extracted MFIs from extract_singlestain_mfi.
#' @param sample_to_rescale data.table of the sample to rescale
#' @param extracted_mfi_df
#' The extracted MFIs from extract_singlestain_mfi. The column "feature" must
#' contain the names of the columns to rescale.
#' By default, if a column of sample_to_rescale is not found in the extracted_mfi_namedlist,
#' the column is rescaled using the minmax method.
#' @param column_negative
#' The column name of the negative population
#' @param column_positive
#' The column name of the positive population
#' @param missing_feature
#' How should missing features (NA) be handled?
#' @param inplace_datatable Should the data.table be modified in place?
#' @param scale_column_fun Function to scale a column
#' @param ...
#' Additional arguments passed to the scale_column_fun
#' @export
rescale_extracted <- function(sample_to_rescale,
                              extracted_mfi,
                              column_negative = c("negative", "unstained"),
                              column_positive = "positive",
                              missing_feature = c("minmax", "center_median"),
                              inplace_datatable = FALSE,
                              scale_column_fun = scale_column_relative, ...) {
    mfis_namedlist <- extract_marker_mfi_list(
        sample_to_rescale = sample_to_rescale,
        extracted_mfi = extracted_mfi,
        column_negative = column_negative,
        column_positive = column_positive
    )
    return(
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
