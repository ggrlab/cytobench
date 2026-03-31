#' Write FCS File to Temporary Memory
#'
#' Attempts to write a `flowFrame` object to a temporary in-memory file.
#' If writing to `/dev/shm/` fails (e.g., on systems without this memory mount),
#' the function falls back to a newly created temporary directory on disk.
#'
#' @param ff A [`flowFrame`](https://rdrr.io/bioc/flowCore/man/flowFrame-class.html) object to be written to a file.
#' @param envir The environment governing the file's lifetime. The temporary file
#'  is deleted when this environment's evaluation completes via [withr::defer()].
#'  Defaults to `rlang::global_env()` for session-long persistence (and
#'  backwards compatibility with the previous `reg.finalizer` approach).
#'  Pass `environment()` to tie the file's lifetime to the calling function's
#'  scope; This is the recommended usage for most workflows:
#'
#'  ```
#'  my_analysis <- function(ff) {
#'      tmp <- write_memory_FCS(ff, envir = environment())
#'      # ... use tmp ...
#'  }
#'  # tmp is cleaned up here
#'  ```
#'
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
write_memory_FCS <- function(ff, envir = rlang::global_env()) {
    ffpath <- base::tempfile()
    if (dir.exists("/dev/shm")) {
        # Use '/dev/shm/{basename(tmpfile)}' instead of the standard tempfile if
        # it exists (shared memory, fast I/O)
        ffpath <- file.path("/dev/shm", basename(ffpath))
    }
    flowCore::write.FCS(ff, ffpath)

    ## Cleanup
    withr::defer(
        unlink(ffpath, recursive = TRUE, force = TRUE),
        envir = envir
    )

    return(ffpath)
}
