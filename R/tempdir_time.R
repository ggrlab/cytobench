#' Create Timestamped Temporary Directory
#'
#' Creates a new temporary directory by appending a timestamp to the default `tempdir()`.
#' This is useful for generating isolated working directories in pipelines or scripts.
#'
#' @return A character string containing the full path to the newly created timestamped directory.
#' The directory is created immediately on disk.
#'
#' @examples
#' tmp <- tempdir_time()
#' list.files(tmp) # should be empty
#'
#' @export
tempdir_time <- function() {
    tmpdir <- tempdir()
    tmpdir_time <- paste0(tmpdir, "_", format(Sys.time(), "%F_%H_%M_%OS2"))
    dir.create(tmpdir_time, recursive = TRUE)
    return(tmpdir_time)
}
