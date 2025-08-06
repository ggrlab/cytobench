#' Export a List of Matrices as CSV Files
#'
#' This function exports a named list of matrices or data frames to `.csv` files,
#' optionally applying a unified feature naming scheme to the column names.
#' Each matrix is saved as a separate CSV file, preserving its name in the output path.
#'
#' @param matrix_list A named list of matrices or data frames to be exported.
#' @param outdir Character. Directory where the CSV files will be saved. Defaults to the current directory (`"."`).
#' @param verbose Logical. Whether to print progress messages during export. Defaults to `TRUE`.
#' @param feature_unified_dict A named data.frame or matrix with rownames corresponding to original
#'   column names and a column `"unified"` providing the standardized feature names.
#'   If `NA`, no renaming is performed.
#'
#' @return Invisibly returns a named character vector of output file paths.
#' @export
#' @keywords cytometry
#' @examples
#' ff <- simulate_ff()
#' export_csv(list("s1.csv" = flowCore::exprs(ff)), outdir = local_tempdir_time())
#'
export_csv <- function(matrix_list,
                       outdir = ".",
                       verbose = TRUE,
                       feature_unified_dict = NA) {
    # Convert each matrix in the list to a data.table for fast writing
    matrix_list_dt <- lapply(matrix_list, data.table::as.data.table)

    # Ensure the output directory exists
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    all_outpaths <- c()

    invisible(lapply(names(matrix_list_dt), function(x) {
        # Create a modified copy with unified column names if a dictionary is provided
        if (!all(is.na(feature_unified_dict))) {
            cp_dt <- data.table::copy(matrix_list_dt[[x]])
            colnames(cp_dt)[colnames(cp_dt) %in% rownames(feature_unified_dict)] <- feature_unified_dict[
                colnames(cp_dt)[colnames(cp_dt) %in% rownames(feature_unified_dict)],
                "unified"
            ]
        } else {
            cp_dt <- matrix_list_dt[[x]]
        }

        # Define output path and create necessary subdirectories
        outpath <- file.path(outdir, x)
        dir.create(dirname(outpath), showWarnings = FALSE, recursive = TRUE)

        # Write CSV file
        data.table::fwrite(cp_dt, outpath)

        # Optionally report progress
        if (verbose) {
            cat("\nWrote ", outpath)
        }

        # Track written file
        all_outpaths <<- c(all_outpaths, outpath)
        names(all_outpaths)[length(all_outpaths)] <<- x
    }))

    if (verbose) {
        cat("\n")
    }

    # Return all output paths invisibly
    invisible(all_outpaths)
}
