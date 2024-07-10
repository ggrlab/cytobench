#' @title Parse filenames into parts
#' @description
#' Extracts defined parts of the filename from the filename string
#' @param fn
#' Filenames, usually a vector of multiple strings
#' @param cn
#' The colnames the filename(s) should be split into
parse_filename <- function(fn,
                           cn = c(
                               "study",
                               "patientID",
                               "tvt",
                               "sample_nr",
                               "processing",
                               "exported_gate",
                               "downsampling_method",
                               "n_cells",
                               "machine",
                               "dye"
                           )) {
    # Extracts our defined parts of the filename
    # from the filename string
    fn_split <- tibble::tibble(
        filename = fn
    ) |>
        dplyr::mutate(
            filename_bn = basename(filename)
        ) |>
        tidyr::separate_wider_delim(filename_bn, delim = "_", names = cn)
    return(
        fn_split
    )
}
