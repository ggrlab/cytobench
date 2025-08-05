#' @title Parse Filenames into Metadata Columns
#'
#' @description
#' Splits structured filenames into multiple metadata components by underscore delimiters.
#' Useful for extracting information such as study ID, sample ID, gating strategy, etc.,
#' from standard naming conventions used in cytometry data processing pipelines.
#'
#' @param fn Character vector. Full paths or base names of filenames to be parsed.
#' @param cn Character vector. Names of the fields to extract from the filename (split by `delim`).
#' Default expects filenames of the format:
#' \preformatted{
#' study_patientID_tvt_sampleNr_processing_exportedGate_downsamplingMethod_nCells_machine_dye.fcs
#' }
#' e.g., \code{"R1_d001_train_s001_raw_cd3_random_1000_navios_KrO.fcs"}.
#' @param delim Character. Delimiter used to split the filename into components.
#' @return A tibble with one row per input filename and columns named by \code{cn},
#' representing each part of the filename, in order. Also includes a column \code{filename}
#' with the original filename.
#'
#' @examples
#' parse_filename(c("R1_d001_train_s001_raw_cd3_random_1000_navios_KrO.fcs"))
#'
#' @export
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
                           ),
                           delim = "_") {
    filename <- filename_bn <- NULL # Linting

    # Convert input to tibble and extract the basename if a path is provided
    fn_split <- tibble::tibble(
        filename = fn
    ) |>
        dplyr::mutate(
            filename_bn = basename(filename)
        ) |>
        tidyr::separate_wider_delim(
            filename_bn,
            delim = delim,
            names = cn
        )

    return(fn_split)
}
