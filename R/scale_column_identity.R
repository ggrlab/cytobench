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
#' ff_dt <- simulate_ff(columns = c("FL1", "FL2", "FL3"), flowcore = FALSE)
#' dt2 <- data.table::data.table(ff_dt)
#' scale_column_identity(
#'     sample_to_rescale = dt2,
#'     scaling_values = c(10, 100),
#'     colX = "FL1",
#'     subtract_bg = TRUE
#' )
#' scale_column_identity(
#'     sample_to_rescale = dt2
#' )
scale_column_identity <- function(sample_to_rescale, ...) {
    invisible(sample_to_rescale)
}
