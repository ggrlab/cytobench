#' @title Update GatingSet Channel and Marker Names
#' @description
#' Helper that renames a `flowWorkspace::GatingSet` using a marker map.
#'
#' @param gs A `flowWorkspace::GatingSet` to rename.
#' @param markermap A `data.frame` with map columns used to rename channels.
#' @param map_oldname
#' A character scalar naming the `markermap` column with current channel names.
#' @param map_newname
#' A character scalar naming the `markermap` column with new channel names.
#' @param map_description
#' A character scalar naming the `markermap` column with marker descriptions.
#'
#' @details
#' The function maps current `GatingSet` channel names to `map_newname` and
#' updates marker names from `map_description`.
#'
#' The map is validated for required columns, complete mappings, and duplicate
#' target channel names. If matching via `map_oldname` fails, it reports the
#' missing channels available mapping names.
#'
#' @return
#' The renamed `flowWorkspace::GatingSet`.
#'
#' @examples
#' \dontrun{
#' gs_new <- gs_update_names(
#'     gs = gs_old,
#'     markermap = marker_map,
#'     map_oldname = "name",
#'     map_newname = "colnames",
#'     map_description = "description_channel.flurochrome.marker.awh"
#' )
#' # Expected: gs_new has renamed channels/markers consistent with marker_map.
#' }
#'
#' @export
gs_update_names <- function(gs, markermap, map_oldname = "name", map_newname = "colnames", map_description = "description_channel.flurochrome.marker.awh") {
    req_cols <- c(map_oldname, map_newname, map_description)
    if (!inherits(gs, "GatingSet")) stop("`gs` must be a flowWorkspace::GatingSet.", call. = FALSE)
    if (!is.data.frame(markermap)) stop("`markermap` must be a data.frame.", call. = FALSE)
    if (!all(req_cols %in% colnames(markermap))) stop("`markermap` is missing required columns.", call. = FALSE)

    old_channel_names <- flowCore::colnames(gs)

    marker_order <- match(old_channel_names, markermap[[map_oldname]])
    if (any(is.na(marker_order))) {
        missing_channels <- old_channel_names[is.na(marker_order)]
        available_names <- unique(as.character(markermap[[map_oldname]]))
        available_names <- available_names[!is.na(available_names)]
        missing_collapsed <- paste(missing_channels, collapse = ", ")
        available_collapsed <- paste(available_names, collapse = ", ")
        message(
            sprintf(
                paste0(
                    "Mapping of old_channel_names failed using `%s`: %d channel(s) not found. \n",
                    "Missing: [%s].\nAvailable: [%s]."
                ),
                map_oldname,
                length(missing_channels),
                missing_collapsed,
                available_collapsed
            )
        )
        stop("Some GatingSet channels are missing from marker map names.", call. = FALSE)
    }

    mapped <- markermap[marker_order, , drop = FALSE]
    new_channel_names <- mapped[[map_newname]]
    new_description_names <- mapped[[map_description]]
    if (any(is.na(new_channel_names)) || anyDuplicated(new_channel_names)) {
        stop("Mapped channel names contain NA or duplicates.", call. = FALSE)
    }
    names(new_description_names) <- new_channel_names

    flowWorkspace::gs_update_channels(
        gs,
        map = data.frame(
            old = old_channel_names,
            new = new_channel_names
        ),
        # Update the flowframe as well?
        # Irrelevant, because we anyways want to apply the gatinghierarchy to new flowframes
        # If `FALSE`, fails.
        all = TRUE
    )
    flowCore::markernames(gs) <- new_description_names
    gs
}

#' @title Check Renamed GatingSet on FlowFrame
#' @description
#' Helper that checks whether a renamed `GatingSet` can be applied to
#' a `flowFrame`, both without and with compensation.
#'
#' @param gs A `flowWorkspace::GatingSet` object.
#' @param ff A `flowCore::flowFrame` used for validation.
#'
#' @return
#' Invisibly returns `TRUE` when both applications succeed.
#'
#' @export
is_updatednames_gs_ff <- function(gs, ff) {
    if (!inherits(gs, "GatingSet")) stop("`gs` must be a flowWorkspace::GatingSet.", call. = FALSE)
    if (!inherits(ff, "flowFrame")) stop("`ff` must be a flowCore::flowFrame.", call. = FALSE)

    tmpfcs <- cytobench::write_memory_FCS(ff[1, ])

    # Validate application without compensation.
    suppressMessages(
        # To suppress 'generating new GatingSet from the gating template...' message
        flowWorkspace::gh_apply_to_cs(
            x = gs,
            cs = flowWorkspace::load_cytoset_from_fcs(tmpfcs),
            compensation_source = "none"
        )
    )

    # Validate application with default compensation handling.
    suppressMessages(flowWorkspace::gh_apply_to_cs(
        x = gs,
        cs = flowWorkspace::load_cytoset_from_fcs(tmpfcs)
        # ,compensation_source = "sample"  ### this is the default!
    ))

    invisible(TRUE)
}
