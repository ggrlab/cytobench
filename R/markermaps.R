#' @title Marker map: Navios DURAClone BASIC
#' @description
#' Static lookup table for mapping channel names and marker labels for the
#' Navios DURAClone BASIC panel.
#' @details
#' This table maps raw channel names (`name`) to harmonized output names
#' (`colnames`) and harmonized marker descriptions
#' (`description_channel.flurochrome.marker.awh`).
#'
#' Intended use is as a `markermap` input to functions that update flowFrame
#' and GatingSet channel/marker names.
#' @format
#' A data.frame with 12 rows and 11 columns:
#' \describe{
#'   \item{panel}{Panel identifier.}
#'   \item{cytometer}{Cytometer identifier.}
#'   \item{name}{Original FCS channel name.}
#'   \item{name_dataset1}{Alternative original name used in dataset1.}
#'   \item{description_short}{Short marker/channel label.}
#'   \item{as_marker}{Logical flag for marker channels.}
#'   \item{fluorochrome}{Fluorochrome label.}
#'   \item{height_width_area}{Signal dimension (A/W).}
#'   \item{channel}{Channel family label.}
#'   \item{colnames}{Harmonized output channel/marker name.}
#'   \item{description_channel.flurochrome.marker.awh}{Combined descriptor.}
#' }
#' @seealso [ff_update_names()], [gs_update_names()]
#' @export
markermap.cytometer_navios.panel_DURACloneBASIC <- data.frame(
    "panel" = "DURACloneBASIC",
    "cytometer" = "Navios",
    "name" = c("FS-A", "FS-W", "SS-A", paste0("FL", c(1:3, 5:8, 10), "-A"), "TIME"),
    "name_dataset1" = c("FS INT", "FS TOF", "SS INT", paste0("FL", c(1:3, 5:8, 10), " INT"), "TIME"),
    "description_short" = c("FSC", "FSC", "SSC", "CD16", "CD56", "CD19", "CD14", "CD4", "CD8", "CD3", "CD45", "TIME"),
    "as_marker" = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
    "fluorochrome" = c(NA, NA, NA, "FITC", "PE", "ECD", "PC7", "APC", "AF700", "AA750", "KrO", NA),
    "height_width_area" = c("A", "W", "A", "A", "A", "A", "A", "A", "A", "A", "A", "NA"),
    "channel" = c(rep("scatter", 3), c(paste0("FL", c(1:3, 5:8, 10))), "time")
) |>
    dplyr::mutate(
        colnames = ifelse(
            is.na(fluorochrome),
            description_short,
            paste0(fluorochrome, "_", description_short)
        ),
        colnames = ifelse(grepl("SC$", colnames), paste0(colnames, "_", height_width_area), colnames),
        description_channel.flurochrome.marker.awh = paste0(
            channel, "-", fluorochrome, "-", description_short,
            ifelse(is.na(height_width_area), "", paste0("-", height_width_area))
        )
    )

#' @title Marker map: Navios DURAClone TCELL
#' @description
#' Static lookup table for mapping channel names and marker labels for the
#' Navios DURAClone TCELL panel.
#' @details
#' This table maps raw channel names (`name`) to harmonized output names
#' (`colnames`) and harmonized marker descriptions
#' (`description_channel.flurochrome.marker.awh`).
#'
#' It also provides `colnames_2025relativisation` for panel-specific downstream
#' relativisation workflows.
#' @format
#' A data.frame with 14 rows and 12 columns:
#' \describe{
#'   \item{panel}{Panel identifier.}
#'   \item{cytometer}{Cytometer identifier.}
#'   \item{name}{Original FCS channel name.}
#'   \item{name_dataset1}{Alternative original name used in dataset1.}
#'   \item{description_short}{Short marker/channel label.}
#'   \item{as_marker}{Logical flag for marker channels.}
#'   \item{fluorochrome}{Fluorochrome label.}
#'   \item{height_width_area}{Signal dimension (A/W).}
#'   \item{channel}{Channel family label.}
#'   \item{colnames}{Harmonized output channel/marker name.}
#'   \item{description_channel.flurochrome.marker.awh}{Combined descriptor.}
#'   \item{colnames_2025relativisation}{Relativisation output label.}
#' }
#' @seealso [ff_update_names()], [gs_update_names()]
#' @export
markermap.cytometer_navios.panel_DURACloneTCELL <- data.frame(
    "panel" = "DURACloneTCELL",
    "cytometer" = "Navios",
    "name" = c("FS-A", "FS-W", "SS-A", paste0("FL", 1:10, "-A"), "TIME"),
    "name_dataset1" = c("FS INT", "FS TOF", "SS INT", paste0("FL", 1:10, " INT"), "TIME"),
    "description_short" = c(
        "FSC", "FSC", "SSC",
        "CD45RA", "CCR7", "CD28", "PD1", "CD27", "CD4", "CD8", "CD3", "CD57", "CD45",
        "TIME"
    ),
    "as_marker" = c(rep(FALSE, 3), rep(TRUE, 10), FALSE),
    "fluorochrome" = c(NA, NA, NA, "FITC", "PE", "ECD", "PC5.5", "PC7", "APC", "AF700", "AA750", "PB", "KrO", NA),
    "height_width_area" = c("A", "W", "A", rep("A", 10), "NA"),
    "channel" = c(rep("scatter", 3), paste0("FL", 1:10), "time")
) |>
    dplyr::mutate(
        colnames = ifelse(
            is.na(fluorochrome),
            description_short,
            paste0(fluorochrome, "_", description_short)
        ),
        colnames = ifelse(grepl("SC$", colnames), paste0(colnames, "_", height_width_area), colnames),
        description_channel.flurochrome.marker.awh = paste0(
            channel, "-", fluorochrome, "-", description_short,
            ifelse(is.na(height_width_area), "", paste0("-", height_width_area))
        ),
        colnames_2025relativisation = paste0(
            ifelse(
                is.na(fluorochrome),
                description_short,
                fluorochrome
            ),
            "-",
            height_width_area
        ),
        colnames_2025relativisation = ifelse(
            grepl("TIME", colnames_2025relativisation),
            "TIME",
            colnames_2025relativisation
        )
    )
