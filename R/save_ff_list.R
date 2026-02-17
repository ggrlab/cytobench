#' Save a Named List of flowFrames as FCS Files
#'
#' Iterates over a named list of flowFrame objects and writes each one
#' to outdir/<name>.fcs, creating directories as needed.
#'
#' @param ff_list A named list of flowFrame objects.
#' @param outdir Directory to write into.
#' @param verbose Logical; print each file path as it is written.
#'
#' @return Invisibly returns the vector of written file paths.
#' @export
#' @examples
#' # Create example flowFrames
#' m <- matrix(runif(100), ncol = 5)
#' colnames(m) <- paste0("Marker", 1:5)
#' ff1 <- flowCore::flowFrame(m)
#' ff2 <- flowCore::flowFrame(m)
#' # Save them to a temporary directory
#' temp_dir <- tempdir()
#' save_ff_list(list(sample1 = ff1, sample2 = ff2), temp_dir)
save_ff_list <- function(ff_list, outdir, verbose = TRUE) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    paths <- character(length(ff_list))
    for (i in seq_along(ff_list)) {
        nm <- names(ff_list)[i]
        outfile <- file.path(outdir, paste0(nm, ".fcs"))
        flowCore::write.FCS(ff_list[[i]], outfile)
        paths[i] <- outfile
        if (isTRUE(verbose)) cat("  Saved:", outfile, "\n")
    }
    invisible(paths)
}
