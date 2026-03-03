#' Rename Channels from Dataset-1 Reference by Position
#'
#' Applies a positional rename strategy where marker names and channel names from
#' a dataset-1 reference frame are copied to a target frame.
#'
#' @param ff A `flowFrame` to rename.
#' @param ds1_reference A reference `flowFrame` (usually from Kaluza dataset 1).
#' @param remove_regex Regex used to normalize reference channel names.
#' @param append_time_marker Logical; append `"TIME"` to marker names.
#'
#' @return Renamed `flowFrame`.
#' @export
rename_by_ds1_position <- function(ff,
                                   ds1_reference,
                                   remove_regex = " ((LIN)|(LOG))",
                                   append_time_marker = TRUE) {
    marker_names <- flowCore::markernames(ds1_reference)
    if (append_time_marker) {
        marker_names <- c(marker_names, "TIME")
    }

    flowCore::markernames(ff) <- structure(unname(marker_names), names = flowCore::colnames(ff))
    flowCore::colnames(ff) <- sub(remove_regex, "", flowCore::colnames(ds1_reference))
    ff
}

#' Build a Channel Name Map Between Two Datasets
#'
#' Creates a data.frame documenting how channel and marker names map between
#' a dataset-1 reference frame, the original dataset-2 frame, and a renamed
#' dataset-2 frame. Useful for auditing rename_by_ds1_position().
#'
#' @param ds1_reference A reference `flowFrame` (e.g. Kaluza dataset 1).
#' @param ds2_original The original (pre-rename) `flowFrame` from dataset 2.
#' @param ds2_renamed The renamed `flowFrame` from dataset 2.
#'
#' @return A `data.frame` with columns `colnames_ds1`, `colnames_ds2`,
#'   `colnames_ds2_original`, `markernames_ds1`, `markernames_ds2`.
#' @export
build_name_map <- function(ds1_reference, ds2_original, ds2_renamed) {
    cbind(
        data.frame(
            colnames_ds1 = flowCore::colnames(ds1_reference),
            colnames_ds2 = flowCore::colnames(ds2_renamed),
            colnames_ds2_original = flowCore::colnames(ds2_original)
        ),
        rbind(
            data.frame(
                markernames_ds1 = flowCore::markernames(ds1_reference),
                markernames_ds2 = flowCore::markernames(ds2_renamed)
            ),
            NA # for TIME marker which has no corresponding marker name
        )
    )
}
