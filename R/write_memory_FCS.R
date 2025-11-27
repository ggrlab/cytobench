#' Write FCS File to Temporary Memory
#'
#' Attempts to write a `flowFrame` object to a temporary in-memory file.
#' If writing to `/dev/shm/` fails (e.g., on systems without this memory mount),
#' the function falls back to a newly created temporary directory on disk.
#'
#' @param ff A [`flowFrame`](https://rdrr.io/bioc/flowCore/man/flowFrame-class.html) object to be written to a file.
#'
#' @return A character string with the file path to the written FCS file.
#'
#' @details
#' This function is useful for workflows requiring in-memory processing or
#' quick temporary storage of FCS files without polluting the working directory.
#'
#' @importFrom flowCore write.FCS
#'
#' @export
#' @keywords cytometry
#' @examples
#' ff <- simulate_ff()
#' tmp <- write_memory_FCS(ff)
#' reread_ff <- flowCore::read.FCS(tmp)
#' identical(reread_ff, ff)
#' nrow(reread_ff) == nrow(ff)
#' print(
#'     summary(flowCore::exprs(ff) - flowCore::exprs(reread_ff))
#' )
write_memory_FCS <- function(ff) {
    tmpfile <- base::tempfile()
    ffpath <- tryCatch(
        {
            # Attempt to write the FCS file to /dev/shm (shared memory, fast I/O)
            ffpath <- file.path("/dev/shm", basename(tmpfile))
            ffpath
        },
        error = function(e) {
            # If /dev/shm is unavailable or writing fails, fallback to a tempfile
            tmpfile
        }
    )

    flowCore::write.FCS(ff, ffpath)
    reg.finalizer(
        e = .GlobalEnv,
        f = function(e) unlink(ffpath, recursive = TRUE, force = TRUE),
        onexit = TRUE
    )
    return(ffpath)
}
