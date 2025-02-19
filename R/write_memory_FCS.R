#' Write FCS File to Memory
#'
#' This function attempts to write a flowFrame object to a temporary file in memory.
#' If writing to the default memory location fails, it creates a new temporary directory
#' and writes the file there.
#'
#' @param ff A flowFrame object to be written to a file.
#' @return A character string representing the path to the written FCS file.
#' @importFrom flowCore write.FCS
#' @examples
#' \dontrun{
#'   library(flowCore)
#'   ff <- read.FCS("path/to/your/file.fcs")
#'   fcs_path <- write_memory.FCS(ff)
#'   print(fcs_path)
#' }
#' @export
write_memory_FCS <- function(ff) {
    tryCatch(
        {
            ffpath <- "/dev/shm/removeme.fcs"
            flowCore::write.FCS(ff, ffpath)
            ffpath
        },
        error = function(e) {
            new_tmpdir <- tempfile()
            dir.create(new_tmpdir)
            ffpath <- file.path(new_tmpdir, "removeme.fcs")
            flowCore::write.FCS(ff, ffpath)
            return(ffpath)
        }
    )
    return(ffpath)
}
