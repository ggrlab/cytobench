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
#' @examples
#' \dontrun{
#' library(flowCore)
#' ff <- read.FCS("path/to/your/file.fcs")
#' fcs_path <- write_memory_FCS(ff)
#' print(fcs_path)
#' }
#'
#' @export
write_memory_FCS <- function(ff) {
    tryCatch(
        {
            # Attempt to write the FCS file to /dev/shm (shared memory, fast I/O)
            ffpath <- "/dev/shm/removeme.fcs"
            flowCore::write.FCS(ff, ffpath)
            return(ffpath)
        },
        error = function(e) {
            # If /dev/shm is unavailable or writing fails, fallback to a temp directory
            new_tmpdir <- tempfile() # Create a unique temporary directory path
            dir.create(new_tmpdir) # Actually create the directory
            ffpath <- file.path(new_tmpdir, "removeme.fcs")
            flowCore::write.FCS(ff, ffpath)
            return(ffpath)
        }
    )
}
