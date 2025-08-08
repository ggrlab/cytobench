#' Create an Empty MFI Tibble
#'
#' Returns a placeholder tibble used for failed or empty MFI extraction.
#' All entries are initialized to `NA`, and column names follow the
#' standard format for single-stain MFI tables.
#'
#' @return A one-row `tibble` with columns:
#'   \describe{
#'     \item{`feature`}{Character column with `NA`.}
#'     \item{`negative`, `positive`}{Numeric columns for median intensities.}
#'     \item{`positive.sd`, `negative.sd`}{Numeric columns for standard deviations.}
#'   }
#' @keywords internal
empty_tibble <- function() {
    tibble::tibble(
        "feature" = NA,
        "negative" = NA,
        "positive" = NA,
        "positive.sd" = NA,
        "negative.sd" = NA
    )
}