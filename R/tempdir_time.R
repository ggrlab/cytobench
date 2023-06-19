#' Create a temporary directory with a timestamp.
#'
#' @export
tempdir_time <- function() {
    tmpdir <- tempdir()
    tmpdir_time <- paste0(tmpdir, "_", format(Sys.time(), "%F_%H_%M_%s"))
    dir.create(tmpdir_time)
    return(tmpdir_time)
}
