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
#' Column name (character) in `df` used for device/sample grouping (default: `"Device"`).
#' @param relevant_columns
#' Character vector of markers to include in the plot. If NULL, all markers are included.
#' @param limit_density_quantile
#' Numeric between 0 and 1; limits plotted densities to a quantile threshold (e.g., 0.95).
#' Useful to suppress extreme density peaks which make all other densities barely visible.
#' Set to `NA` to disable.
#'
#' @return A `ggplot2` object showing density curves across markers, devices, and samples.
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

    # Use all columns (channels/markers) if no subset is given
    if (all(is.null(relevant_columns))) {
        relevant_columns <- flowCore::colnames(ff_gated[[1]])
    }

    # Convert single transform function to named list if needed
    transformlist <- transformlist_named(transformlist, relevant_columns)
    # Apply identity function to unused columns to keep them unchanged
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

    # Merge metadata: File, Device, Sample
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
#' Additional arguments passed to [calc_densities()] if `densities_dt` is not already in
#' density format (i.e. missing `x` and `y` columns). This allows users to pass grouping
#' and transformation parameters directly when providing raw expression data. See
#'  [calc_densities()] for details on accepted arguments.
#'
#' @return A [ggplot2::ggplot] object.
#' @examples
#' library(data.table)
#'
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
