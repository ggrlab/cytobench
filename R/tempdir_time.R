#' Create Timestamped Temporary Directory
#'
#' Creates a new temporary directory by appending a timestamp to the default `tempdir()`.
#' This is useful for generating isolated working directories in pipelines or scripts.
#' @param dir
#' Base directory to append the timestamp to. Defaults to [tempdir()]. You might want
#' to use `withr::local_tempdir()` to create a temporary directory that is automatically
#' cleaned up after the session ends (Useful for testing).
#' @param ... Additional arguments passed to `dir.create()`, such as `showWarnings` or `recursive`.
#' @return A character string containing the full path to the newly created timestamped directory.
#' The directory is created immediately on disk.
#'
#' @examples
#' tmp <- tempdir_time()
#' list.files(tmp) # should be empty
#'
#' @export
tempdir_time <- function(dir = tempdir(), ...) {
    tmpdir_time <- paste0(dir, "_", format(Sys.time(), "%F_%H_%M_%OS3"))
    dir.create(tmpdir_time, ...)
    return(tmpdir_time)
}

#' Create Timestamped Subdirectory Inside a Local Temporary Directory
#'
#' Generates a timestamped subdirectory within a temporary directory that will be
#' automatically cleaned up when the current environment exits (e.g., during testing).
#' This is useful for managing temporary file outputs in isolated, reproducible sessions.
#' @param showWarnings Logical. If `TRUE`, shows warnings when creating the directory.
#' From `base::dir.create()`.
#' @param ... Additional arguments passed to `tempdir_time()`.
#'
#' @return A character string containing the path to the created timestamped subdirectory.
#'
#' @details
#' Uses [withr::local_tempdir()] to create a base temporary directory that is scoped to the
#' current environment (typically cleaned up on exit), and then calls `tempdir_time()` to
#' create a unique subdirectory within it.
#'
#' @examples
#' \dontrun{
#' tmp <- local_tempdir_time()
#' print(tmp)
#' }
#'
#' @export
local_tempdir_time <- function(showWarnings = FALSE, ...) {
    # Create an environment-scoped temporary base directory
    tmpdir <- withr::local_tempdir(.local_envir = parent.frame(2))

    # Create a timestamped subdirectory inside the base temp dir
    tempdir_time(file.path(tmpdir, "dirAutoremove"), showWarnings = showWarnings, ...)
}
