#' Plot Density Distributions for Gated Flow Cytometry Data
#'
#' This function generates density plots for multiple flow cytometry markers across different devices and samples.
#' It applies optional transformations, merges metadata, and creates a faceted plot of density curves.
#'
#' Every (unique!) flowframe must have a unique name present in `df$File`, which will be
#' matched to the `"Sample"` column. One sample would often (not necessarily!)
#' correspond to multiple flowframes, e.g., if multiple files were acquired
#' for the same sample on multiple devices.
#'
#' @param ff_gated
#' Named list of `flowFrame` objects, each representing gated flow cytometry data for one sample.
#' Names must correspond to entries in `df$File`. Can be a flowset or a named list of flowFrames.
#' @param df
#' A `data.table` containing metadata for samples. Must include `"File"`, `dfcol_grouping_samples`, and `"Sample"` columns.
#' @param device_colors
#' Optional named vector of colors for each device (e.g., `c("Fortessa" = "blue", "Aurora" = "red")`).
#' @param transformlist Either a single transformation function or a named list of functions, one per marker.
#'                      If NULL, no transformation is applied.
#' @param density_n Integer; number of points for kernel density estimation (default: 500).
#' @param dfcol_grouping_samples
#' Character scalar naming the column in `df` used for grouping/coloring
#' samples (default: `"Device"`).
#' @param relevant_columns
#' Character vector of markers to include in the plot. If NULL, all markers are included.
#' @param limit_density_quantile
#' Numeric between 0 and 1; limits plotted densities to a quantile threshold (e.g., 0.95).
#' Useful to suppress extreme density peaks which make all other densities barely visible.
#' Set to `NA` to disable.
#'
#' @return A `ggplot2` object showing density curves across markers, devices, and samples.
#'
#' @details
#' The function expects that every element of `ff_gated` has a unique name and that
#' these names match entries in `df$File`. Densities are computed per `File`, then
#' joined with metadata (`File`, grouping column, `Sample`) for faceting and coloring.
#' @export
#' @examples
#' fs <- cytobench::simulate_fs(n_samples = 3, ncells = 250, columns = c("CD4", "CD8"))
#'
#' # Metadata
#' df_meta <- data.table::data.table(
#'     File = flowCore::sampleNames(fs),
#'     Device = c("A", "B", "A"),
#'     Sample = c("S1", "S2", "S2")
#' )
#'
#' # Colors for devices
#' device_cols <- c("A" = "steelblue", "B" = "firebrick")
#'
#' # Apply density plot function
#' plot <- plot_densities_ff(
#'     ff_gated = fs,
#'     df = df_meta,
#'     device_colors = device_cols,
#'     transformlist = function(x) asinh(x / 1000),
#'     relevant_columns = c("CD4", "CD8")
#' )
#'
#' # Show plot
#' print(plot)
plot_densities_ff <- function(
  ff_gated,
  df,
  device_colors,
  transformlist = NULL,
  density_n = 500,
  dfcol_grouping_samples = "Device",
  relevant_columns = NULL,
  limit_density_quantile = NA
) {
    . <- variable <- value <- File <- y <- x <- V1 <- NULL # R CMD check compatibility
    if ("flowSet" %in% class(ff_gated)) {
        ff_gated <- flowCore::flowSet_to_list(ff_gated)
    }
    # Ensure ff_gated is a named list
    if (!all(!is.null(names(ff_gated)))) {
        stop("ff_gated must be a named list of flowFrame objects")
    }
    if (all(is.null(relevant_columns))) {
        relevant_columns <- flowCore::colnames(ff_gated[[1]])
    }

    # Backward-compatible second guard: use all columns when no marker subset is given
    if (all(is.null(relevant_columns))) {
        relevant_columns <- flowCore::colnames(ff_gated[[1]])
    }

    # Normalize transform input to a named list keyed by marker
    transformlist <- cyCompare::transformlist_named(transformlist, relevant_columns)
    # Keep non-selected columns unchanged via identity transforms
    for (col_x in flowCore::colnames(ff_gated[[1]])[!flowCore::colnames(ff_gated[[1]]) %in% relevant_columns]) {
        transformlist[[col_x]] <- identity
    }


    # Extract expression data and melt into long format
    gated_dt <- lapply(ff_gated, function(x) {
        flowCore::exprs(x) |> data.table::as.data.table()
    }) |>
        data.table::rbindlist(idcol = "File", fill = TRUE)

    densities <- calc_densities(
        dt = gated_dt,
        groupings = "File",
        markers = relevant_columns,
        transformlist = transformlist,
        limit_density_quantile = limit_density_quantile
    )

    # Keep only metadata fields needed for annotation in the final plot
    df_part <- data.table::data.table(df)
    df_part <- df_part[, c("File", dfcol_grouping_samples[[1]], "Sample"), with = FALSE]
    densities <- densities[df_part, on = "File"]


    # Plot densities by Sample and Marker (x = intensity, y = density)
    p0 <- plot_densities(
        densities_dt = densities,
        column_marker = "marker",
        column_color = dfcol_grouping_samples[[1]]
    ) +
        ggh4x::facet_grid2(Sample ~ marker, scales = "free", independent = "y")
    # Add manual colors if provided
    if (!all(is.null(device_colors))) {
        p0 <- p0 +
            ggplot2::scale_color_manual(values = device_colors) +
            ggplot2::scale_fill_manual(values = device_colors)
    }

    return(p0)
}
#' Compute a 1D kernel density estimate
#'
#' Computes a kernel density estimate for a numeric vector, optionally applying
#' a transformation prior to density estimation. Missing values are removed
#' silently. If no finite values remain, a single-row NA result is returned.
#'
#' @param x Numeric vector of expression values.
#' @param transform_x Optional transformation function applied to `x` prior to
#'   density estimation (e.g. `asinh`). Defaults to the identity.
#' @param density_n Integer specifying the number of points used for the density
#'   estimate (passed to [stats::density()]).
#'
#' @return A [data.table::data.table] with columns:
#' \describe{
#'   \item{x}{Evaluation points of the density.}
#'   \item{y}{Estimated density values.}
#' }
#'
#' @examples
#' x <- rnorm(1000)
#' d <- compute_density(x)
#' head(d)
#'
#' # With transformation
#' d_asinh <- compute_density(x, transform_x = asinh)
#' head(d_asinh)
#'
#' @export
compute_density <- function(x, transform_x = NULL, density_n = 500) {
    x_noNA <- x[!is.na(x)]

    if (length(x_noNA) == 0) {
        return(data.table::data.table(x = NA_real_, y = NA_real_))
    }

    if (is.null(transform_x)) {
        transform_x <- identity
    }

    d <- stats::density(transform_x(x_noNA), n = density_n)

    data.table::data.table(
        x = d$x,
        y = d$y
    )
}


#' Calculate expression density estimates per group and marker
#'
#' Converts a wide expression table into long format and computes kernel density
#' estimates for each combination of grouping variables and markers.
#'
#' @param dt A data.frame or data.table containing expression values.
#' @param groupings Character vector of column names used for grouping
#'   (e.g. `"File"`). May be `NULL` for ungrouped densities.
#' @param markers Character vector of marker columns. Defaults to all columns
#'   except `groupings`.
#' @param transformlist Optional named list of transformation functions, indexed
#'   by marker name from `markers`. Missing entries default to no transformation.
#' @param limit_density_quantile Optional numeric in `(0, 1]`. If provided,
#'   density values are capped at this quantile within each marker to suppress
#'   extreme peaks.
#' @param density_n Integer specifying the number of points used for each kernel
#'   density estimate (passed to [compute_density()]).
#'
#' @return A [data.table::data.table] with columns:
#' \describe{
#'   \item{groupings}{Grouping variables (if provided).}
#'   \item{marker}{Marker name.}
#'   \item{x}{Evaluation points of the density.}
#'   \item{y}{Estimated density values.}
#' }
#'
#' @examples
#' library(data.table)
#'
#' # Toy expression table
#' dt <- data.table(
#'     File = rep(c("A", "B"), each = 100),
#'     CD4  = rnorm(200),
#'     CD8  = rnorm(200, mean = 1)
#' )
#'
#' dens <- calc_densities(
#'     dt,
#'     groupings = "File"
#' )
#'
#' dens[marker == "CD4"][1:5]
#'
#' @export
calc_densities <- function(dt,
                           groupings = NULL,
                           markers = colnames(dt)[!colnames(dt) %in% groupings],
                           transformlist = NULL,
                           limit_density_quantile = NA,
                           density_n = 500) {
    dt_long <- data.table::melt(
        dt[, c(groupings, markers), with = FALSE],
        id.vars = groupings,
        measure.vars = markers,
        variable.name = "marker",
        value.name = "expression"
    )

    if (is.character(dt_long[["expression"]])) {
        stop(
            "Expression values after melting cannot be character. ",
            "Did you set groupings and markers correctly?"
        )
    }

    densities <- dt_long[
        ,
        {
            compute_density(
                expression,
                transformlist[[as.character(marker[[1]])]],
                density_n = density_n
            )
        },
        by = c(groupings, "marker")
    ]

    if (!is.na(limit_density_quantile)) {
        top_quantile <- densities[
            ,
            stats::quantile(y, limit_density_quantile),
            by = "marker"
        ]

        densities <- densities[
            top_quantile,
            on = "marker"
        ][
            ,
            y := pmin(y, V1)
        ][
            ,
            V1 := NULL
        ]
    }
    densities
}

#' Plot expression density curves
#'
#' Visualizes kernel density estimates produced by [calc_densities()] using
#' faceting by marker and optional grouping by color.
#'
#' @param densities_dt A data.table as returned by [calc_densities()].
#' @param column_marker Character scalar giving the column used for faceting
#'   markers (default: `"marker"`).
#' @param column_color Optional character scalar giving a column used for
#'   coloring and grouping densities (e.g. `"File"`).
#' @param ...
#' Additional arguments passed to [calc_densities()] when `densities_dt` is not
#' already in density format (i.e. missing `x` and `y` columns). This allows
#' users to pass grouping and transformation parameters directly when providing
#' raw expression data. See [calc_densities()] for accepted arguments.
#'
#' @return A ggplot2::ggplot object.
#' If x/y are missing, interpret input as raw expression data and compute densities.
#' @examples
#' dt <- data.table(
#'     File = rep(c("A", "B"), each = 200),
#'     CD4  = rnorm(400),
#'     CD8  = rnorm(400, mean = 0.5)
#' )
#'
#' dens <- calc_densities(dt, groupings = "File")
#'
#' # Colored by File
#' plot_densities(dens, column_color = "File")
#'
#' # Uncolored (aggregated view). This looks super weird because the densities
#' # are computed per File, but then all plotted together without grouping.
#' plot_densities(dens)
#'
#' @export
plot_densities <- function(densities_dt,
                           column_marker = "marker",
                           column_color = NULL,
                           ...) {
    if (!all(c("x", "y") %in% colnames(densities_dt))) {
        # Then I ASSUME that calc_densities has not been called yet, and that the user
        # is passing the original long data.table with expression values. So I will call
        # calc_densities here.
        densities_dt <- calc_densities(
            dt = densities_dt,
            groupings = column_color,
            ...
        )
    }
    if (is.null(column_color)) {
        aes_ggplot <- ggplot2::aes(x = x, y = y)
        aes_ribbon <- ggplot2::aes(ymin = 0, ymax = y)
        rows <- NULL
    } else {
        aes_ggplot <- ggplot2::aes(
            x = x,
            y = y,
            color = !!rlang::sym(column_color)
        )
        aes_ribbon <- ggplot2::aes(
            ymin = 0,
            ymax = y,
            fill = !!rlang::sym(column_color)
        )
        rows <- ggplot2::vars(!!rlang::sym(column_color))
    }
    remaining_cols <- setdiff(
        colnames(densities_dt),
        c("x", "y", column_marker, column_color)
    )
    if (length(remaining_cols) > 0) {
        warning(
            "The following columns in densities_dt are not used for plotting and may lead to wrong/weird plots: ",
            paste(remaining_cols, collapse = ", ")
        )
    }

    ggplot2::ggplot(densities_dt, aes_ggplot) +
        ggplot2::geom_line() +
        ggplot2::geom_ribbon(aes_ribbon, alpha = 0.2) +
        ggh4x::facet_grid2(
            cols = ggplot2::vars(!!rlang::sym(column_marker)),
            rows = rows,
            scales = "free",
            independent = "y"
        ) +
        ggpubr::theme_pubr() +
        ggplot2::theme(
            axis.line.y  = ggplot2::element_blank(),
            axis.text.y  = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
        ) +
        ggplot2::ylab("Density") +
        ggplot2::xlab("Transformed MFI")
}



#' Plot expression density curves
#'
#' Visualizes kernel density estimates produced by [calc_densities()] as
#' multiple base R plots, split by marker and optional color group.
#'
#' @param densities_dt A data.table as returned by [calc_densities()].
#' @param column_marker Character scalar giving the column used to split
#'   markers into separate plots (default: `"marker"`).
#' @param column_color Optional character scalar giving a column used for
#'   coloring and grouping densities (e.g. `"File"`).
#' @param color_map Optional named character vector of colors used when
#'   `column_color` is provided.
#' @param ribbon_alpha Numeric in `[0, 1]` controlling ribbon transparency.
#' @param main_prefix Optional character scalar prepended to each plot title.
#' @param ...
#' Additional arguments passed to [calc_densities()] when `densities_dt` is not
#' already in density format (i.e. missing `x` and `y` columns). This allows
#' users to pass grouping and transformation parameters directly when providing
#' raw expression data. See [calc_densities()] for accepted arguments.
#'
#' @return Invisibly returns a list with metadata about the generated plots.
#' If x/y are missing, interpret input as raw expression data and compute densities.
#' @examples
#' dt <- data.table(
#'     File = rep(c("A", "B"), each = 200),
#'     CD4  = rnorm(400),
#'     CD8  = rnorm(400, mean = 0.5)
#' )
#'
#' dens <- calc_densities(dt, groupings = "File")
#'
#' # Colored by File
#' plot_densities_base(dens, column_color = "File")
#'
#' # Uncolored (aggregated view). This looks super weird because the densities
#' # are computed per File, but then all plotted together without grouping.
#' plot_densities_base(dens)
#'
#' @export
plot_densities_base <- function(densities_dt,
                                column_marker = "marker",
                                column_color = NULL,
                                color_map = NULL,
                                ribbon_alpha = 0.2,
                                main_prefix = NULL,
                                ...) {
    if (!all(c("x", "y") %in% colnames(densities_dt))) {
        # Then I ASSUME that calc_densities has not been called yet, and that the user
        # is passing the original long data.table with expression values. So I will call
        # calc_densities here.
        densities_dt <- calc_densities(
            dt = densities_dt,
            groupings = column_color,
            ...
        )
    }

    if (!column_marker %in% colnames(densities_dt)) {
        stop("column_marker not found in densities_dt: ", column_marker)
    }
    if (!is.null(column_color) && !column_color %in% colnames(densities_dt)) {
        stop("column_color not found in densities_dt: ", column_color)
    }


    if (is.null(column_color)) {
        color_map <- c("all" = "black")
    } else {
        color_values <- unique(as.character(densities_dt[[column_color]]))
        if (is.null(color_map)) {
            color_map <- stats::setNames(
                grDevices::hcl.colors(length(color_values), palette = "Dark 3"),
                color_values
            )
        }
    }

    marker_values <- unique(as.character(densities_dt[[column_marker]]))
    panel_count <- 0L
    for (marker_x in marker_values) {
        for (color_x in names(color_map)) {
            if (color_x == "all") {
                # All the same color
                dt_plot <- densities_dt[get(column_marker) == marker_x]
                panel_title <- paste0(column_marker, ": ", marker_x)
                current_color <- "black"
            } else {
                dt_plot <- densities_dt[get(column_marker) == marker_x & get(column_color) == color_x]
                panel_title <- paste0(column_marker, ": ", marker_x, " | ", column_color, ": ", color_x)
                current_color <- color_map[[color_x]]
            }

            if (!is.null(main_prefix)) {
                panel_title <- paste0(main_prefix, " | ", panel_title)
            }

            plotted <- plot_density_base_single(
                dt_plot = dt_plot,
                panel_title = panel_title,
                color = current_color,
                ribbon_alpha = ribbon_alpha
            )
            if (isTRUE(plotted)) {
                panel_count <- panel_count + 1L
            }
        }
    }

    invisible(
        list(
            n_plots = panel_count,
            markers = marker_values,
            column_color = column_color
        )
    )
}


#' Draw a single base-R density panel
#'
#' Helper used by [plot_densities_base()] to render one marker/group
#' combination. The function validates finite coordinates, initializes an empty
#' plotting area, draws a filled ribbon, and overlays the density line.
#'
#' @param dt_plot A `data.frame`-like object with numeric `x` and `y` columns.
#' @param panel_title Character scalar used as plot title.
#' @param color Color used for the density line and ribbon fill.
#' @param ribbon_alpha Numeric in `[0, 1]` controlling ribbon transparency.
#'
#' @return Logical scalar indicating whether a panel was drawn successfully.
#' @export
plot_density_base_single <- function(dt_plot, panel_title = NULL, color = "black", ribbon_alpha = .2) {
    dt_plot <- dt_plot[is.finite(dt_plot[["x"]]) & is.finite(dt_plot[["y"]]), , drop = FALSE]
    if (nrow(dt_plot) == 0) {
        warning("Skipping panel with no finite x/y values: ", panel_title)
        return(FALSE)
    }
    dt_plot <- dt_plot[order(dt_plot[["x"]]), , drop = FALSE]

    x_vals <- dt_plot[["x"]]
    y_vals <- dt_plot[["y"]]
    y_max <- max(y_vals, na.rm = TRUE)
    if (!is.finite(y_max) || y_max <= 0) {
        y_max <- 1
    }

    graphics::plot(
        x_vals,
        y_vals,
        type = "n",
        ylim = c(0, y_max),
        xlab = "Transformed MFI",
        ylab = "Density",
        main = panel_title,
        yaxs = "i",
        yaxt = "n"
    )
    graphics::polygon(
        x = c(x_vals, rev(x_vals)),
        y = c(rep(0, length(x_vals)), rev(y_vals)),
        col = grDevices::adjustcolor(color, alpha.f = ribbon_alpha),
        border = NA
    )
    graphics::lines(x_vals, y_vals, col = color, lwd = 2)
    TRUE
}
