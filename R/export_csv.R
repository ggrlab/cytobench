#' @export 
export_csv <- function(
    matrix_list,
    outdir = ".",
    verbose = TRUE,
    feature_unified_dict = NA) {
    matrix_list_dt <- lapply(matrix_list, data.table::as.data.table)
    dir.create(outdir, recursive = TRUE)
    invisible(lapply(names(matrix_list), function(x) {
        # The following lines are only such that the .csv files have
        # always the same unified column names, regardless if it
        # was empty (--> single staining sample) or not
        #
        # But i do not want to change the data.table itself, so i copy it
        # after export_fcs writes the "original" column names into the "description"
        if (!all(is.na(feature_unified_dict))) {
            cp_dt <- data.table::copy(matrix_list[[x]])
            colnames(cp_dt)[colnames(cp_dt) %in% rownames(f_dict)] <- f_dict[
                colnames(cp_dt)[colnames(cp_dt) %in% rownames(f_dict)],
                "unified"
            ]
        } else {
            cp_dt <- matrix_list[[x]]
        }
        outpath <- file.path(outdir, x)
        dir.create(dirname(outpath), showWarnings = FALSE, recursive = TRUE)
        data.table::fwrite(cp_dt, outpath)
        if (verbose) {
            cat("\nWrote ", outpath)
        }
    }))
    if (verbose) {
        cat("\n")
    }
}
