#' Rescale a sample using a named list of extracted MFIs
#'
#' This function rescales all columns of a sample using a named list of extracted MFIs.
#' @param sample_to_rescale data.table of the sample to rescale
#' @param extracted_mfi_namedlist
#'    The extracted MFIs as list of lists:
#         $`FITC-A`
#         $`FITC-A`$negative
#         [1] 495.8738
#         $`FITC-A`$positive
#         function (v)
#         .approxfun(x, y, v, method, yleft, yright, f, na.rm)
#         <bytecode: 0x55e899638948>
#         <environment: 0x55e89f0ca338>
#'    In contrast to rescale_named(), this is a named list of LISTs with the names "positive" and
#'    "negative" MFIs or have no name.
#'    names(extracted_mfi_namedlist) should be the names of the columns in sample_to_rescale.
#'    names(extracted_mfi_namedlist[1]) should either be "positive", "negative" or have no name
#'    If it has no name, I assume that this column should be regarded as a "missing feature".
#'
#'    By default, if a column of sample_to_rescale is not found in the extracted_mfi_namedlist,
#'    the column is rescaled using the minmax method.
#' @param missing_feature
#' How should missing features (NA) be handled?
#' @param inplace_datatable Should the data.table be modified in place?
#' @param scale_column_fun Function to scale a column
#' @param ...
#' Additional arguments passed to the scale_column_fun
#' @export
rescale_named_quantiles <- function(sample_to_rescale,
                                    extracted_mfi_namedlist,
                                    reference_marker,
                                    missing_feature = c("minmax", "center_median", NA),
                                    inplace_datatable = FALSE,
                                    scale_column_fun = scale_column_minmaxQ, ...) {
    # tmp <- extracted_mfi_namedlist
    # extracted_mfi_namedlist <- tmp
    if (!inplace_datatable || !data.table::is.data.table(sample_to_rescale)) {
        # That here effectively generates a new datatable
        # If the matrix _IS_ already a datatable, effectively a copy is created.
        # If the matrix is NOT created, this is _always_ called
        sample_to_rescale <- data.table::as.data.table(sample_to_rescale)
    }

    missing_mfis <- colnames(sample_to_rescale)[!colnames(sample_to_rescale) %in% names(extracted_mfi_namedlist)]
    if (length(missing_mfis) > 0) {
        missing_mfis_list <- rep(list(list(NA)), length(missing_mfis))
        names(missing_mfis_list) <- missing_mfis
        extracted_mfi_namedlist <- c(extracted_mfi_namedlist, missing_mfis_list)
    }
    if (length(missing_mfis) > 0 && !is.na(missing_feature[1])) {
        warning(paste0(
            "No MFI feature values for \'", missing_mfis, "\', added NA to the extracted_mfi_namedlist\n"
        ))
    }

    given_scale_column_fun <- scale_column_fun
    # ECDF: Empirical Cumulative Distribution Function
    sample_to_rescale[, reference_quantiles := stats::ecdf(get(reference_marker))(get(reference_marker))]
    # data.table::setkeyv(sample_to_rescale, "reference_quantiles")

    for (colX in colnames(sample_to_rescale)) {
        scale_column_fun <- given_scale_column_fun
        if (all(is.na(extracted_mfi_namedlist[[colX]]))) {
            if (is.na(missing_feature[1])) {
                # Then do not rescale this column
                extracted_values <- NA
            } else if (missing_feature[1] == "center_median") {
                extracted_values <- sample_to_rescale[, lapply(.SD, median), .SDcols = colX]
                scale_column_fun <- scale_column_minmax
            } else if (missing_feature[1] == "minmax") {
                extracted_values <- c(
                    sample_to_rescale[, lapply(.SD, max), .SDcols = colX],
                    sample_to_rescale[, lapply(.SD, min), .SDcols = colX]
                )
                scale_column_fun <- scale_column_minmax
            } else {
                stop("Missing '", colX, "' rescaling")
            }
        } else {
            extracted_values <- extracted_mfi_namedlist[[colX]]
        }
        if (all(is.na(extracted_values))) {
            # Then do not rescale this column
            next
        }
        extracted_values <- unlist(extracted_values)

        scale_column_fun(
            sample_to_rescale = sample_to_rescale,
            scaling_values = extracted_values,
            colX = colX,
            ...
        )
    }
    return(sample_to_rescale)
}
