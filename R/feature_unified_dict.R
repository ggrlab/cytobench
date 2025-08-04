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
