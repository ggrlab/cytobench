#' Apply Marker-Specific Transformation to a List of Cytometry Data Tables
#'
#' Applies an asinh or user-defined transformation to each column of each data table in a list,
#' based on provided marker-specific scaling factors. Intended for flow cytometry data pre-processing.
#'
#' @param datatable_list A list of `data.table` objects containing raw cytometry data.
#' @param scalings_names_unified
#' A named numeric vector assigning a scaling factor to each marker name (e.g., `"FL1 CD45RA FITC" = 50000`).
#'   If a column has value `notransform_special_number`, no transformation is applied to it.
#' @param notransform_special_number
#' Numeric sentinel (default: -14) indicating that no transformation should be applied to the column.
#' @param feature_unified_dict
#' A data.frame with a row for each feature and a `"unified"` column used to map alternative column names
#'   to the canonical names found in `scalings_names_unified`. Typically obtained via `feature_unified_dict_default()`.
#' @param transform_fun
#' A function with signature `function(x_values, scaling_factor)` used to apply the transformation.
#'   Default is `asinh(x / scaling_factor)`.
#' @param inplace Logical. If `TRUE`, the function modifies the input data.tables in-place. If `FALSE`, a copy is made.
#'
#' @return The transformed list of data.tables. Columns not listed in `scalings_names_unified` will raise an error.
#'
#' @details
#' If `inplace = FALSE`, copies of the input data.tables are created before transformation. If `inplace = TRUE`, the
#' transformations are done in-place using `data.table::set()`.
#' The function will fail if unknown columns or NAs are encountered.
#'
#' @export
#' @keywords cytometry
#' @examples
# Simulate a list of data.tables with some flow-like markers
#' dt1 <- data.table::data.table(
#'     `FL1 CD45RA FITC` = runif(100, 0, 100000),
#'     `FL2 CCR7 PE` = runif(100, 0, 100000),
#'     `FS INT` = runif(100, 0, 500)
#' )
#'
#' dt2 <- data.table::data.table(
#'     `FL1 CD45RA FITC` = runif(100, 0, 100000),
#'     `FL2 CCR7 PE` = runif(100, 0, 100000),
#'     `FS INT` = runif(100, 0, 500)
#' )
#'
#' datatable_list <- list(dt1, dt2)
#'
#' # Apply transformation
#' transformed_list <- transform_dts(lapply(datatable_list, data.table::data.table))
#'
#' # Inspect the result
#' summary(transformed_list[[1]][, `FL1 CD45RA FITC`])
#' summary(datatable_list[[1]][, `FL1 CD45RA FITC`])
transform_dts <- function(datatable_list,
                          scalings_names_unified = c(
                              "FS INT" = -14,
                              "FS TOF" = -14,
                              "SS INT" = -14,
                              "FL1 CD45RA FITC" = 50e3,
                              "FL2 CCR7 PE" = 50e3,
                              "FL3 CD28 ECD" = 50e3,
                              "FL4 PD1 PC5.5" = 50e3,
                              "FL5 CD27 PC7" = 50e3,
                              "FL6 CD4 APC" = 50e3,
                              "FL7 CD8 AF700" = 50e3,
                              "FL8 CD3 AA750" = 50e3,
                              "FL9 CD57 PB" = 50e3,
                              "FL10 CD45 KrO" = 50e3,
                              "TIME" = -14
                          ),
                          notransform_special_number = -14,
                          feature_unified_dict = feature_unified_dict_default(),
                          transform_fun = function(x_values, scaling_factor) {
                              asinh(x_values / scaling_factor)
                          },
                          inplace = TRUE) {
    # If not modifying in-place, convert each element to data.table
    if (!inplace) {
        datatable_list <- lapply(datatable_list, data.table::as.data.table)
    }

    # Apply transformation to each data.table in the list
    lapply(datatable_list, function(cytoX) {
        # Ensure there are no NA values
        if (any(is.na(cytoX))) {
            stop("Why are there NAs in the sample?")
        }

        # Loop through each column
        for (col_i in seq_len(ncol(cytoX))) {
            col_name <- colnames(cytoX)[col_i]

            # Determine the unified name for the current column
            if (col_name %in% rownames(feature_unified_dict)) {
                ss_unified <- feature_unified_dict[col_name, ][["unified"]]
                col_scaling <- scalings_names_unified[ss_unified]
            } else {
                col_scaling <- scalings_names_unified[col_name]
            }

            # Handle transformation decision
            if (is.na(col_scaling)) {
                stop(paste0(col_name, ": Unknown scaling"))
            } else if (col_scaling == notransform_special_number) {
                # No transformation for this marker
                next
            } else {
                # Apply transformation
                data.table::set(
                    cytoX,
                    j = col_i,
                    value = transform_fun(cytoX[, col_i, with = FALSE], col_scaling)
                )
            }
        }
    })

    return(datatable_list)
}
