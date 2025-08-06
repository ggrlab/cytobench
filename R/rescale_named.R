#' Rescale a sample using a named list of MFIs
#'
#' This function rescales each column in a sample data table using a named list of extracted marker intensities (MFIs).
#' Each column is matched to its corresponding MFI vector (e.g., with named elements `"negative"` and `"positive"`).
#'
#' @param sample_to_rescale A `data.table` representing a single sample. Each column corresponds to a marker/channel.
#'
#' @param extracted_mfi_namedlist
#' A named list where each name corresponds to a column in `sample_to_rescale`.
#' Each element is a numeric vector representing the MFI values for that marker/channel.
#' If a column is missing in this list, it will be handled according to `missing_feature`.
#' See also extracted_mfi in `extract_marker_mfi_list()`.
#'
#' @param missing_feature Character. How to handle missing MFI entries (i.e., where list entry is `NA`).
#' Options are:
#' \itemize{
#'   \item `"minmax"` (default): Rescale using column-wise min/max.
#'   \item `"center_median"`: Shift column so that median is centered.
#'   \item `NA`: Do not rescale columns with missing MFI values.
#' }
#'
#' @param inplace_datatable Logical. If `TRUE`, modifies `sample_to_rescale` in place. If `FALSE`, returns a new copy. Default is `FALSE`.
#'
#' @param scale_column_fun Function. The scaling function to apply per column. Default is `scale_column_minmax()`.
#' This function must accept at least `sample_to_rescale` (data.table), `scaling_values` (vector of numbers), and `colX` (character, a column of `sample_to_rescale`).
#' @param ... Additional arguments passed to `scale_column_fun`.
#'
#' @return A `data.table` with rescaled columns. If `inplace_datatable = TRUE`, the same object is modified and returned invisibly.
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
#' # Convert extracted MFI table into named list of feature -> MFI (negative/positive)
#' mfis_namedlist <- extract_marker_mfi_list(
#'     sample_to_rescale = flowCore::exprs(fs_ss[["sample0_12-panel"]]),
#'     extracted_mfi = extracted_mfis_singlestain,
#'     column_negative = "negative",
#'     column_positive = "positive"
#' )
#'
#' rescale_named(
#'     sample_to_rescale = flowCore::exprs(fs_ss[["sample0_12-panel"]]),
#'     extracted_mfi_namedlist = mfis_namedlist,
#'     missing_feature = "minmax",
#'     inplace_datatable = TRUE,
#'     scale_column_fun = scale_column_minmax
#' )
rescale_named <- function(sample_to_rescale,
                          extracted_mfi_namedlist,
                          missing_feature = c("minmax", "center_median"),
                          inplace_datatable = FALSE,
                          scale_column_fun = scale_column_minmax,
                          ...) {
    # Ensure we work with a data.table (copy if needed)
    if (!inplace_datatable || !data.table::is.data.table(sample_to_rescale)) {
        # That here effectively generates a new datatable
        # If the matrix _IS_ already a datatable, effectively a copy is created.
        # If the matrix is NOT created, this is _always_ called
        sample_to_rescale <- data.table::as.data.table(sample_to_rescale)
    }

    # Identify columns missing in the MFI list and append NA entries
    missing_mfis <- colnames(sample_to_rescale)[!colnames(sample_to_rescale) %in% names(extracted_mfi_namedlist)]
    if (length(missing_mfis) > 0) {
        missing_mfis_list <- rep(list(NA), length(missing_mfis))
        names(missing_mfis_list) <- missing_mfis
        extracted_mfi_namedlist <- c(extracted_mfi_namedlist, missing_mfis_list)

        if (length(missing_mfis) > 0 && !is.na(missing_feature[1])) {
            warning(paste0(
                "No MFI feature values for: ", paste(missing_mfis, collapse = ", "),
                "; added NA to extracted_mfi_namedlist."
            ))
        }
    }

    # Save original scaling function to allow fallback override
    given_scale_column_fun <- scale_column_fun

    # Iterate over each column to apply scaling
    for (colX in colnames(sample_to_rescale)) {
        scale_column_fun <- given_scale_column_fun

        # Determine how to obtain or approximate scaling values
        if (all(is.na(extracted_mfi_namedlist[[colX]]))) {
            if (is.na(missing_feature[1])) {
                # Then do not rescale this column
                extracted_values <- NA
            } else if (missing_feature[1] == "center_median") {
                extracted_values <- sample_to_rescale[, lapply(.SD, stats::median), .SDcols = colX]
                scale_column_fun <- scale_column_minmax
            } else if (missing_feature[1] == "minmax") {
                extracted_values <- c(
                    sample_to_rescale[, lapply(.SD, max), .SDcols = colX],
                    sample_to_rescale[, lapply(.SD, min), .SDcols = colX]
                )
                scale_column_fun <- scale_column_minmax
            } else {
                stop("Unknown missing_feature strategy for column: ", colX)
            }
        } else {
            extracted_values <- extracted_mfi_namedlist[[colX]]
        }

        if (all(is.na(extracted_values))) {
            # Then do not rescale this column
            next
        }

        # Convert to unlisted numeric vector if needed
        extracted_values <- unlist(extracted_values)

        # Apply scaling function
        scale_column_fun(
            sample_to_rescale = sample_to_rescale,
            scaling_values = extracted_values,
            colX = colX,
            ...
        )
    }

    return(sample_to_rescale)
}
