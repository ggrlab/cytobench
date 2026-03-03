#' @title Update `flowFrame` Channel and Marker Names
#' @description
#' Internal helper that renames a `flowFrame` using a mapping table and keeps
#' related compensation/autofluorescence keywords synchronized.
#'
#' @param flowframe A `flowCore::flowFrame` to rename.
#' @param markermap A `data.frame` with columns `name`, `colnames`, and
#' `description_channel.flurochrome.marker.awh`.
#'
#' @details
#' The function matches current `flowFrame` channel names against
#' `markermap$name`, updates channel names (`colnames`) and marker names
#' (`markernames`), then rewrites spillover/autofluorescence keyword names.
#' It validates required map columns, missing mappings, and duplicate targets.
#'
#' @return
#' A renamed `flowCore::flowFrame` with updated channel, marker, and keyword
#' metadata.
#'
#' @examples
#' \dontrun{
#' ff_new <- ff_update_names(flowframe = ff_old, markermap = marker_map)
#' # Expected: ff_new has renamed channels/markers and aligned spillover names.
#' }
#'
#' @keywords internal
ff_update_names <- function(flowframe,
                            markermap = markermap.cytometer_navios.panel_DURACloneTCELL[, c(
                                map_oldname, map_newname, map_description
                            )],
                            spillmat_names = NULL,
                            map_oldname = "name",
                            map_newname = "colnames",
                            map_description = "description_channel.flurochrome.marker.awh") {
    req_cols <- c(map_oldname, map_newname, map_description)
    if (!inherits(flowframe, "flowFrame")) {
        stop("`flowframe` must be a flowCore::flowFrame.", call. = FALSE)
    }
    if (!is.data.frame(markermap)) {
        stop("`markermap` must be a data.frame.", call. = FALSE)
    }
    if (!all(req_cols %in% colnames(markermap))) {
        stop("`markermap` is missing required columns.", call. = FALSE)
    }

    old_channel_names <- flowCore::colnames(flowframe)
    marker_order <- match(old_channel_names, markermap[[map_oldname]])
    if (any(is.na(marker_order))) {
        stop("Some flowFrame channels are missing from `markermap$name`.",
            call. = FALSE
        )
    }

    markermap <- tibble::as_tibble(markermap)
    mapped <- markermap[marker_order, ]
    new_channel_names <- mapped[[map_newname]]
    new_marker_names <- mapped[[map_description]]
    if (any(is.na(new_channel_names)) || anyDuplicated(new_channel_names)) {
        stop("Mapped channel names contain NA or duplicates.", call. = FALSE)
    }
    names(new_marker_names) <- new_channel_names

    ff_old <- flowframe[1, ]
    flowCore::colnames(flowframe) <- new_channel_names
    flowCore::markernames(flowframe) <- new_marker_names

    # Update spillover matrices in keywords and validate compensation syntax.
    if (is.null(spillmat_names)) {
        spillmats <- tryCatch(
            flowCore::spillover(flowframe),
            error = function(e) {
                warning("Could not extract spillover matrices: ", e$message)
                list()
            }
        )
        spillmat_names <- names(spillmats)
    }
    spillmats <- flowCore::keyword(flowframe)[spillmat_names]

    for (spillmat_name in names(spillmats)) {
        spillmat <- spillmats[[spillmat_name]]
        if (is.null(spillmat)) next
        new_parameternames <- get_new_parameternames(
            ff_old, spillmat,
            new_channel_names, markermap[[map_oldname]]
        )
        colnames(spillmats[[spillmat_name]]) <- new_parameternames
        rownames(spillmats[[spillmat_name]]) <- NULL

        # test the compensation:
        flowCore::compensate(flowframe[1, ], spillover = spillmats[[spillmat_name]])
        # Replace the spillover matrix in the flowframe keywords
        flowCore::keyword(flowframe)[[spillmat_name]] <- spillmats[[spillmat_name]]
    }

    # Update autofluorescence keyword vectors to renamed channel names.
    keyword_names <- names(flowCore::keyword(flowframe))
    kw_af <- grep("autofluorescence", keyword_names, ignore.case = TRUE, value = TRUE)
    for (kw_af_single in kw_af) {
        af <- flowCore::keyword(flowframe)[[kw_af_single]]
        names(af) <- new_parameternames
        flowCore::keyword(flowframe)[[kw_af_single]] <- af
    }

    # Remove malformed/internal helper keys from keywords.
    kw_names <- names(flowCore::keyword(flowframe))
    flowCore::keyword(flowframe)[is.na(kw_names)] <- NULL
    flowCore::keyword(flowframe)[grepl("flowCore", kw_names)] <- NULL
    flowframe
}

#' @title Map Spillover Parameter Names to Renamed Channels
#' @description
#' Internal helper that maps spillover matrix parameter names from the original
#' `flowFrame` channel naming to the renamed channel naming.
#'
#' @param ff A `flowCore::flowFrame` containing original `$PnN` keywords.
#' @param spillmat A spillover matrix with either "numeric" (numbers as character) column
#' names (channel indices) or original channel-name column names.
#' @param new_channel_names A character vector of renamed channel names aligned
#' to `markermap_names`.
#' @param markermap_names A character vector of original channel names used as
#' matching reference.
#'
#' @details
#' The function reads old parameter names from `ff` using `$P{index}N` keyword
#' lookups when `colnames(spillmat)` are numeric indices. If `spillmat`
#' already has channel-name columns, those names are used directly. In both
#' cases, names are mapped to `new_channel_names` via matching against
#' `markermap_names`.
#'
#' To check if the colnames are numeric indices, if ANY colnames result in NA
#' when coerced to numeric, we assume that ALL of them are NOT numeric indices but
#' should be treated literally.
#'
#' @return
#' A character vector of renamed parameter names, in the same order as
#' `colnames(spillmat)`.
#'
#' @examples
#' \dontrun{
#' new_pn <- get_new_parameternames(
#'     ff = ff,
#'     spillmat = spillmat,
#'     new_channel_names = c("CD3", "CD4"),
#'     markermap_names = c("FL1", "FL2")
#' )
#' # Expected: a character vector aligned with spillmat parameter order.
#' }
#'
#' @export
#'
get_new_parameternames <- function(ff, spillmat, new_channel_names, markermap_names) {
    # Example spillmats:
    #     $`$SPILLOVER`
    #           4     5
    #  [1,] 1.000 0.191
    #  [2,] 0.008 1.000
    # $`$SPILLOVER`
    #      FITC_CD16 PE_CD56
    # [1,]     1.000   0.186
    # [2,]     0.009   1.000

    # Get the old parameternames to map the spillover matrices
    spill_cols <- colnames(spillmat)
    if (is.null(spill_cols)) {
        stop("`spillmat` must have column names.", call. = FALSE)
    }
    is_indexed <- all(grepl("^[0-9]+$", spill_cols))
    if (is_indexed) {
        old_parameternames <- flowCore::keyword(ff)[paste0("$P", spill_cols, "N")] |> unlist()
    } else {
        old_parameternames <- spill_cols
    }
    new_parameternames <- new_channel_names[match(old_parameternames, markermap_names)]

    if (length(new_parameternames) != ncol(spillmat)) {
        stop("Length of new parameter names does not match number of spillover columns.")
    }
    if (any(is.na(new_parameternames))) {
        stop("Could not map one or more spillover parameter names.", call. = FALSE)
    }
    new_parameternames
}
