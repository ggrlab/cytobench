#' Create Default Unified Feature Mapping Table
#'
#' This function returns a data frame used to harmonize marker/channel names across different
#' flow cytometry naming conventions. It defines a set of canonical (unified) feature names
#' and known variations (alternatives), useful for unifying data across single-stain, panel,
#' and empty/placeholder conditions.
#'
#' The returned data frame can be used to rename columns in expression matrices or metadata
#' tables, e.g., for CSV or FCS exports.
#'
#' @return A `data.frame` with the following columns:
#' \describe{
#'   \item{`unified`}{Canonical marker name (e.g., `"FL1 CD45RA FITC"`).}
#'   \item{`unified_single_staining`}{Simplified name used for single-staining formats (e.g., `"FL1 CD3 FITC"`).}
#'   \item{`alt_i`}{Identifier for which alternative format the row represents (`"alt0"`, `"alt1"`, `"alt2"`).}
#'   \item{`alternative`}{A variant of the unified name (e.g., `"FL1 CD3 FITC"`, `"FL1 empty FITC"`).}
#' }
#'
#' Row names of the returned object are set to the `alternative` values for easy lookup.
#'
#' @examples
#' dict <- feature_unified_dict_default()
#' dict["FL1 CD3 FITC", "unified"]
#' @export
feature_unified_dict_default <- function() {
    alternative <- unified <- NULL # lint
    f_dict <- data.frame(unified = c(
        # "FS INT", "SS INT", "FS TOF", "TIME",
        "FL1 CD45RA FITC",
        "FL2 CCR7 PE",
        "FL3 CD28 ECD",
        "FL4 PD1 PC5.5",
        "FL5 CD27 PC7",
        "FL6 CD4 APC",
        "FL7 CD8 AF700",
        "FL8 CD3 AA750",
        "FL9 CD57 PB",
        "FL10 CD45 KrO"
    )) |>
        dplyr::mutate(
            alt0 = unified,
            alt1 = sub("(FL[0-9]+) [^ ]* ", "\\1 CD3 ", unified),
            unified_single_staining = sub("(FL[0-9]+) [^ ]* ", "\\1 CD3 ", unified),
            alt2 = sub("(FL[0-9]+) [^ ]* ", "\\1 empty ", unified)
        ) |>
        tidyr::pivot_longer(
            cols = tidyr::starts_with("alt"),
            names_to = "alt_i", values_to = "alternative"
        ) |>
        dplyr::distinct(alternative, .keep_all = TRUE) |>
        as.data.frame()

    rownames(f_dict) <- f_dict$alternative
    return(f_dict)
}
