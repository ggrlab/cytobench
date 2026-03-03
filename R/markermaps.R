# For each cytometer panel, it is necessary to have a marker map that maps:
#  - `name`: The column names in the FCS files
# to
#  - `colnames`:
#       The new column names in the output data
#  - `description_channel.flurochrome.marker.awh`:
#       The new description including channel, fluorochrome, marker, and area/width/height split by "-"
#  - `as_marker`: Whether this channel should be regarded as a marker or not (e.g., FSC, SSC, Time are not markers)

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
