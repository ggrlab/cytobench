#' Base R engine for pairwise marker scatterplots
#'
#' This function creates pairwise scatterplots of markers in a data frame using base R plotting functions.
#'
#' @inheritParams plot_markers_pairwise
#' @param col
#' Color for points in scatterplots (default: semi-transparent black). Is overridden when using `geom = "pointdensity"`.
#' @param parargs
#' `engine = "base"`. A list of graphical parameters to pass to `par()` for base R plotting.
#' @param mgp
#' `engine = "base"`. A numeric vector of length 3 specifying the margin line for the axis title, labels, and line.
#' @param method.args
#' A list of additional arguments to pass to the density calculation method in `ggpointdensity:::StatPointdensity$compute_group()`.
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
                                           mfrow = c(ncol(df_mat), ncol(df_mat)) - 1,
                                           mar = c(3, 3, 0, 0),
                                           oma = c(0, 0, 0, 0),
                                           xaxs = "i", yaxs = "i"
                                       ),
                                       mgp = c(2, 1, 0),
                                       geom = c("points", "pointdensity"),
                                       col = grDevices::rgb(0, 0, 0, .2),
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

    markernames <- colnames(df_mat)
    par(parargs)
    # Loop through each pair and produce a plot
    for (marker_y in markernames[-1]) {
        for (marker_x in markernames) {
            if (verbose) {
                cat(marker_x, "vs", marker_y, ": ")
            }
            xy_str <- paste0(marker_x, "_", marker_y)
            if (marker_x == marker_y) {
                if (verbose) {
                    cat("diagonal\n")
                }
                # next
            } else if (!xy_str %in% all_combos_str) {
                if (verbose) {
                    cat("empty\n")
                }
                plot.new()
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
            }
        }
    }
}
