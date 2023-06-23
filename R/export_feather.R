#' @export
export_feather <- function(matrix_list,
                           outdir = ".",
                           verbose = TRUE,
                           feature_unified_dict = NA,
                           add_feather_ending = TRUE) {
    matrix_list_dt <- lapply(matrix_list, data.table::as.data.table)
    dir.create(outdir, recursive = TRUE)
    invisible(lapply(names(matrix_list_dt), function(x) {
        # The following lines are only such that the .csv files have
        # always the same unified column names, regardless if it
        # was empty (--> single staining sample) or not
        #
        # But i do not want to change the data.table itself, so i copy it
        # after export_fcs writes the "original" column names into the "description"
        if (!all(is.na(feature_unified_dict))) {
            cp_dt <- data.table::copy(matrix_list_dt[[x]])
            colnames(cp_dt)[colnames(cp_dt) %in% rownames(feature_unified_dict)] <- feature_unified_dict[
                colnames(cp_dt)[colnames(cp_dt) %in% rownames(feature_unified_dict)],
                "unified"
            ]
        } else {
            cp_dt <- matrix_list_dt[[x]]
        }
        outpath <- file.path(outdir, x)
        dir.create(dirname(outpath), showWarnings = FALSE, recursive = TRUE)
        if (add_feather_ending) {
            outpath <- paste0(outpath, ".feather")
        }
        feather::write_feather(cp_dt, outpath)
        if (verbose) {
            cat("\nWrote ", outpath)
        }
    }))
    if (verbose) {
        cat("\n")
    }
}
