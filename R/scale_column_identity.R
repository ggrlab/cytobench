#' @title Identity Scaling Function (No Operation)
#'
#' @description
#' This function performs no scaling on the data. It is used as a placeholder or for debugging
#' when no transformation is required. It invisibly returns the input unchanged.
#'
#' @param sample_to_rescale A `data.table` (or compatible object) to be returned unchanged.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns `sample_to_rescale` without modification.
#'
#' @export
#' @keywords relativisation
scale_column_identity <- function(sample_to_rescale, ...) {
    invisible(sample_to_rescale)
}
