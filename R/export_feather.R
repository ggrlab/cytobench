#' Export a List of Matrices as Feather Files
#'
#' This function exports a named list of matrices or data frames to the Apache Feather file format.
#' Optionally, it can rename columns according to a unified feature dictionary before export.
#'
#' @param matrix_list A named list of matrices or data frames to be exported.
#' @param outdir Character. Path to the output directory. Defaults to the current directory.
#' @param verbose Logical. Whether to print messages during export. Defaults to TRUE.
#' @param feature_unified_dict A data.frame or matrix with rownames corresponding to original feature names
#'   and a column named `"unified"` containing the unified feature names. Used to rename columns before export.
#'   Set to `NA` to disable renaming.
#' @param add_feather_ending Logical. Whether to append ".feather" to output filenames. Defaults to TRUE.
#'
#' @return Invisibly returns a named character vector of output file paths.
#' @export
export_feather <- function(matrix_list,
                           outdir = ".",
                           verbose = TRUE,
                           feature_unified_dict = NA,
                           add_feather_ending = TRUE) {
    
    # Convert all list elements to data.tables
    matrix_list_dt <- lapply(matrix_list, data.table::as.data.table)
    
    # Ensure the output directory exists
    dir.create(outdir, recursive = TRUE)
    
    all_outpaths <- c()
    
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
        
        # Construct output path
        outpath <- file.path(outdir, x)
        dir.create(dirname(outpath), showWarnings = FALSE, recursive = TRUE)
        
        # Add ".feather" file ending if specified
        if (add_feather_ending) {
            outpath <- paste0(outpath, ".feather")
        }
        
        # Write data to feather format
        feather::write_feather(cp_dt, outpath)
        
        # Optionally print status
        if (verbose) {
            cat("\nWrote ", outpath)
        }
        
        # Save the output path in a named vector
        all_outpaths <<- c(all_outpaths, outpath)
        names(all_outpaths)[length(all_outpaths)] <<- x
    }))
    
    if (verbose) {
        cat("\n")
    }
    
    invisible(all_outpaths)
}
