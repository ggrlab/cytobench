#' @export
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
                          }, inplace = TRUE) {
    if (!inplace) {
        # That here effectively generates a new datatable
        # If the matrix _IS_ already a datatable, effectively a copy is created.
        # If the matrix is NOT created, this is _always_ called
        matrix_list <- lapply(matrix_list, data.table::as.data.table)
    }
    only_null_results <- lapply(matrix_list, function(cytoX) {
        if (any(is.na(cytoX))) {
            stop("Why are there NAs in the sample?")
        }
        for (col_i in seq_len(ncol(cytoX))) {
            if (colnames(cytoX)[col_i] %in% rownames(feature_unified_dict)) {
                ss_unified <- feature_unified_dict[colnames(cytoX)[col_i], ][["unified"]]
                col_scaling <- scalings_names_unified[ss_unified]
            } else {
                col_scaling <- scalings_names_unified[colnames(cytoX)[col_i]]
            }

            if (is.na(col_scaling)) {
                stop(paste0(colnames(cytoX)[col_i], ": Unknown scaling"))
            } else if (col_scaling == notransform_special_number) {
                # do nothing
            } else {
                data.table::set(
                    cytoX,
                    j = col_i,
                    value = transform_fun(cytoX[, col_i, with = FALSE], col_scaling)
                )
            }
        }
    })
    return(matrix_list)
}
