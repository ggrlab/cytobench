#' Base R engine for pairwise marker scatterplots
#'
#' This function creates pairwise scatterplots of markers in a data frame using base R plotting functions.
#'
#' @inheritParams plot_markers_pairwise
#' @param col
#' Color for points in scatterplots (default: semi-transparent black). Is overridden when using `geom = "pointdensity"`.
#' @param parargs
#' `engine = "base"`. A list of graphical parameters to pass to `par()` for base R plotting in the
#' INDIVIDUAL scatterplots! Take care not to override layout parameters - this will likely lead to
#' wrongly arranged plots.
#' @param mgp
#' `engine = "base"`. A numeric vector of length 3 specifying the margin line for the axis title, labels, and line.
#' @param firstcol_width
#' `engine = "base"`. Width of the first column (marker names) relative to other columns (scatterplots).
#' @param firstrow_height
#' `engine = "base"`. Height of the first row (marker names) relative to other rows (scatterplots).
#' @param cex_title
#' `engine = "base"`. Character expansion factor for the global title.
#' @param method.args
#' A list of additional arguments to pass to the density calculation method in
#' `ggpointdensity:::StatPointdensity$compute_group()`.
#' See `?ggpointdensity::geom_pointdensity` for details.
#' @param adjust
#' See `?ggpointdensity::geom_pointdensity` for details.
#' @param n_colors
#' Number of colors to use for density coloring in pointdensity geom (default: 256).
#' @param ... Additional arguments passed to the plotting function `scattermore::scattermoreplot()`.
plot_markers_pairwise_base <- function(df,
                                       cofactor_namedvec,
                                       special_cofactor_list = list(),
                                       transform_fun = asinh,
                                       transform_fun_name = "asinh",
                                       verbose = FALSE,
                                       parargs = list(
                                           mar = c(3, 3, 0, 0),
                                           xaxs = "i", yaxs = "i"
                                       ),
                                       mgp = c(2, 1, 0),
                                       geom = c("points", "pointdensity"),
                                       col = grDevices::rgb(0, 0, 0, .2),
                                       title_global = NULL,
                                       firstcol_width = 0.15,
                                       firstrow_height = 0.15,
                                       cex_title = 2,
                                       modelines = TRUE,
                                       # arguments for pointdensity
                                       method.args = list(),
                                       adjust = 1,
                                       count_transform = function(x) x,
                                       n_colors = 256,
                                       ...) {
    if (geom[1] == "hex") {
        geom <- "points"
        warning("geom = 'hex' not implemented in base engine; using 'points' instead")
    }
    # All pairwise combinations of markers
    all_combos <- utils::combn(colnames(df), 2, simplify = FALSE)
    all_combos_str <- lapply(all_combos, paste0, collapse = "_")

    df_mat <- as.matrix(df)
    if (!missing(cofactor_namedvec)) {
        df_mat[, names(cofactor_namedvec)] <- df_mat[, names(cofactor_namedvec)] %*% (diag(1 / cofactor_namedvec))
        df_mat <- transform_fun(df_mat)
    }
    if (modelines) {
        if (!require(modeest)) {
            stop("Package 'modeest' is required for adding mode lines. Please install it or set modelines = FALSE.")
        }
        modes <- apply(df_mat, 2, modeest::mlv1)
    }

    markernames <- colnames(df_mat)
    layoutmat <- matrix(0, nrow = length(markernames) + 1, ncol = length(markernames) + 1)
    rownames(layoutmat) <- c("namecol", markernames)
    colnames(layoutmat) <- c("namecol", markernames)
    n <- colnames(layoutmat)
    plotcount <- 0
    for (j in 0:length(markernames)) {
        for (i in 0:length(markernames)) {
            if (i == j || i == 1) {
                next
            }
            if ((n[i + 1] == "namecol" || n[j + 1] == "namecol")) {
                plotcount <- plotcount + 1
                layoutmat[i + 1, j + 1] <- plotcount
            } else if (j < i) {
                plotcount <- plotcount + 1
                layoutmat[i + 1, j + 1] <- plotcount
            }
        }
    }
    layoutmat <- layoutmat[-2, ]
    layoutmat <- layoutmat[, -ncol(layoutmat)]

    # > layoutmat
    #         namecol FITC-A PE-A
    # namecol       0      3    6
    # PE-A          1      4    0
    # ECD-A         2      5    7

    layout(mat = layoutmat, widths = c(firstcol_width, rep(1, length(markernames) - 1)), heights = c(firstrow_height, rep(1, length(markernames) - 1)))
    # layout.show(max(layoutmat))

    plot_n <- 0
    for (marker_x in colnames(layoutmat)) {
        for (marker_y in rownames(layoutmat)) {
            par(parargs)
            if (layoutmat[marker_y, marker_x] <= 0) {
                next
            }
            if (marker_x == "namecol") {
                plot_n <- plot_n + 1
                # First column: marker names
                par(mar = c(0, 0, 0, 0))
                plot.new()
                text(0.5, 0.5, marker_y, cex = 2, srt = 90)
            } else if (marker_y == "namecol") {
                plot_n <- plot_n + 1
                par(mar = c(0, 0, 0, 0))
                plot.new()
                text(0.5, 0.5, marker_x, cex = 2)
            } else {
                if (verbose) {
                    cat(marker_x, "vs", marker_y, ": ")
                }
                xy_str <- paste0(marker_x, "_", marker_y)
                if (marker_x == marker_y || !xy_str %in% all_combos_str) {
                } else {
                    if (verbose) {
                        cat("plotting\n")
                    }
                    if (xy_str %in% names(special_cofactor_list)) {
                        vals_x <- transform_fun(df[[marker_x]] / special_cofactor_list[[xy_str]][1])
                        vals_y <- transform_fun(df[[marker_y]] / special_cofactor_list[[xy_str]][2])

                        xylab <- lapply(
                            c("x" = marker_x, "y" = marker_y),
                            function(marker) {
                                paste0(
                                    transform_fun_name, "(", marker, "/",
                                    special_cofactor_list[[xy_str]][[ifelse(marker == marker_x, 1, 2)]], "))"
                                )
                            }
                        )
                    } else {
                        vals_x <- df_mat[, marker_x]
                        vals_y <- df_mat[, marker_y]
                        xylab <- paste0(
                            transform_fun_name, "(", c(marker_x, marker_y), "/",
                            cofactor_namedvec[c(marker_x, marker_y)], ")"
                        )
                        names(xylab) <- c("x", "y")
                    }
                    plot_n <- plot_n + 1
                    if (geom[1] == "points") {
                        scattermore::scattermoreplot(
                            vals_x, vals_y,
                            xlab = xylab["x"],
                            ylab = xylab["y"],
                            frame.plot = TRUE,
                            mgp = mgp,
                            col = col,
                            ...
                        )
                    } else if (geom[1] == "pointdensity") {
                        if (length(vals_x) > 20000) {
                            method_densitycalc <- "kde2d"
                        } else {
                            method_densitycalc <- "neighbors"
                        }
                        points_density <- ggpointdensity:::StatPointdensity$compute_group(
                            data.frame(x = vals_x, y = vals_y),
                            # compute_group() uses get_limits() and dimension() of ggplot scales.
                            # Therefore, implement dummy scales here:
                            scales = list(
                                x = list(get_limits = function() range(vals_x), dimension = function() range(vals_x)),
                                y = list(get_limits = function() range(vals_y), dimension = function() range(vals_y))
                            ),
                            method = method_densitycalc,
                            method.args = method.args,
                            adjust = adjust
                        )
                        points_density_transformed <- count_transform(points_density$density)
                        colfunc <- grDevices::colorRampPalette(viridis::viridis(2))(n_colors)
                        # https://github.com/LKremer/ggpointdensity/blob/master/R/geom_pointdensity.R
                        scattermore::scattermoreplot(
                            vals_x, vals_y,
                            xlab = xylab["x"],
                            ylab = xylab["y"],
                            frame.plot = TRUE,
                            mgp = mgp,
                            # col = colfunc(count_transform(points_density$density)),
                            col = colfunc[cut(points_density_transformed, breaks = n_colors)],
                            ...
                        )
                    } else {
                        stop("Unknown geom: ", geom[1])
                    }
                    if (modelines) {
                        abline(v = modes[marker_x], h = modes[marker_y], col = "red")
                    }
                }
            }
        }
    }
    mtext(title_global, side = 3, line = -5, outer = TRUE, cex = cex_title)
}
