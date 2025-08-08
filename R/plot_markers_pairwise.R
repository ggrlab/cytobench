#' Plot all pairwise marker combinations with given cofactors
#'
#' Generates a matrix of 2D-density plots for all marker combinations in a flow cytometry dataset.
#' Axes are asinh-transformed and scaled using marker-specific cofactors. This is particularly useful
#' for visual QC of compensation, cofactor choices and signal separation.
#'
#' @param ff A `flowFrame` containing the cytometry data.
#' @param cofactor_namedvec
#' A named numeric vector providing default cofactors per marker. The names
#' of the vector should be the same as the marker names in the flowFrame object.
#' @param special_cofactor_list
#' A named list with overrides for specific marker combinations.
#' Each name should be `"markerX_markerY"`, with values being a
#' numeric vector of length 2: `c(cofactor_x, cofactor_y)`.
#' @param transform_fun A function used to transform intensity values (default: `asinh`).
#' @param transform_fun_name Character label for the transformation function (used in axis labels).
#' @param geom Character.
#' One of `"hex"` (default), `"points"`, or `"pointdensity"` -
#' determines the ggplot2 geometry used. Hex is usually still the fastest.
#' @param bins Number of bins (resolution) in each 2D plot (default: 50).
#' @param diag_plot Logical. Whether to include labeled marker names along the diagonal (default: FALSE).
#' @param debugplots Logical. If TRUE, plots text placeholders instead of actual data - useful for layout testing.
#' @param axis_full_labels Logical. If TRUE, axis labels include marker names and cofactor formulas.
#' @param n_cells Number of cells to downsample from `ff` (default: `Inf` = use all).
#' @param count_transform Function to transform hexbin/point counts (e.g. `log10(count + 1)`).
#' @param add_ggplot_elements A list of ggplot layers to add to each panel (e.g. `list(ggplot2::theme_minimal())`).
#'
#' @return A `patchwork` object: a composite grid of all marker pairwise comparisons.
#'
#' @export
#'
#' @examples
#' fs <- simulate_fs(
#'     n_samples = 2,
#'     flowcore = TRUE,
#'     ncells = 250,
#'     columns = c("FL1-A", "FL2-A", "FL3-A")
#' )
#' p0 <- plot_markers_pairwise(fs[[1]])
#' print(p0)
#' p1 <- plot_markers_pairwise(
#'     fs[[1]],
#'     cofactor_namedvec = c("FL1-A" = 5, "FL2-A" = 10, "FL3-A" = 15),
#'     special_cofactor_list = list("FL1-A_FL2-A" = c(3, 6), "FL2-A_FL3-A" = c(4, 8)),
#'     transform_fun = function(x) log10(x + 1),
#'     transform_fun_name = "log10p1"
#' )
#' print(p1)
#' p2 <- plot_markers_pairwise(
#'     fs[[1]],
#'     cofactor_namedvec = c("FL1-A" = 5, "FL2-A" = 10, "FL3-A" = 15),
#'     special_cofactor_list = list("FL1-A_FL2-A" = c(3, 6), "FL2-A_FL3-A" = c(4, 8)),
#'     geom = "points",
#'     bins = 100,
#'     diag_plot = TRUE,
#'     debugplots = FALSE,
#'     axis_full_labels = TRUE,
#'     n_cells = 1000
#' )
#' print(p2)
#' p3 <- plot_markers_pairwise(
#'     fs[[1]],
#'     cofactor_namedvec = c("FL1-A" = 5, "FL2-A" = 10, "FL3-A" = 15),
#'     special_cofactor_list = list("FL1-A_FL2-A" = c(3, 6), "FL2-A_FL3-A" = c(4, 8)),
#'     geom = "pointdensity",
#'     bins = 100,
#'     diag_plot = TRUE,
#'     debugplots = FALSE,
#'     axis_full_labels = TRUE,
#'     n_cells = 1000,
#'     add_ggplot_elements = list(ggplot2::theme_minimal())
#' )
#' print(p3)
#' @keywords cytometry
plot_markers_pairwise <- function(ff,
                                  cofactor_namedvec,
                                  special_cofactor_list,
                                  transform_fun = asinh,
                                  transform_fun_name = "asinh",
                                  geom = c("hex", "points", "pointdensity"),
                                  bins = 50,
                                  diag_plot = FALSE,
                                  debugplots = FALSE,
                                  axis_full_labels = TRUE,
                                  n_cells = Inf,
                                  count_transform = function(x) log10(x + 1),
                                  add_ggplot_elements = list()) {
    x <- y <- density <- count <- NULL # LINTR
    # Generate blank diagonal marker name plots
    marker_names <- names(flowCore::markernames(ff))
    xy_plots <- list()
    xy_plots_rotated <- list()
    for (marker_x in marker_names) {
        xy_plots[[marker_x]] <- ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0, y = 0, size = 18, label = marker_x) +
            ggplot2::theme_void()
        xy_plots_rotated[[marker_x]] <- ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0, y = 0, size = 18, label = marker_x, angle = 90) +
            ggplot2::theme_void()
    }

    # All pairwise combinations of markers
    all_combos <- utils::combn(marker_names, 2, simplify = FALSE)
    all_combos_str <- lapply(all_combos, paste0, collapse = "_")

    # Extract data and optionally downsample
    exprs_dt <- data.table::as.data.table(flowCore::exprs(ff))
    if (!is.infinite(n_cells)) {
        exprs_dt <- exprs_dt[sample(.N, min(.N, n_cells))]
    }

    # Start assembling plot matrix
    plots_all_patchwork <- list()
    for (marker_y in marker_names) {
        plots_all_patchwork <- c(plots_all_patchwork, xy_plots_rotated[marker_y])
        for (marker_x in marker_names) {
            xy_str <- paste0(marker_x, "_", marker_y)

            which_matching <- which(unlist(all_combos_str) == xy_str)
            # Only plot lower triangle of matrix
            if (xy_str %in% all_combos_str) {
                if (!missing(cofactor_namedvec)) {
                    cofactor_x <- cofactor_namedvec[[marker_x]]
                    cofactor_y <- cofactor_namedvec[[marker_y]]
                } else {
                    cofactor_x <- 1
                    cofactor_y <- 1
                }
                if (!missing(special_cofactor_list) && xy_str %in% names(special_cofactor_list)) {
                    cofactor_x <- special_cofactor_list[[xy_str]][1]
                    cofactor_y <- special_cofactor_list[[xy_str]][2]
                }

                if (debugplots) {
                    p <- ggplot2::ggplot() +
                        ggplot2::annotate("text",
                            x = 0, y = 0, size = 15,
                            label = paste0(xy_str, "\nx: ", cofactor_x, " y: ", cofactor_y)
                        ) +
                        ggplot2::theme_void()
                } else {
                    dt <- tibble::tibble(
                        x = transform_fun(exprs_dt[[marker_x]] / cofactor_x),
                        y = transform_fun(exprs_dt[[marker_y]] / cofactor_y)
                    )
                    p <- ggplot2::ggplot(dt, ggplot2::aes(x = x, y = y))

                    if (geom[1] == "hex") {
                        irrelevant <- hexbin::BTC(2) # this is a completely irrelevant piece of code
                        # but geom_hex depends on the hexbin package.
                        # If I "import" hexbin, it tells me that it's not used (and therefore a NOTE)
                        # So, I here use it in an irrelevant way to avoid the NOTE while
                        # still importing it
                        p <- p + ggplot2::geom_hex(
                            bins = bins,
                            ggplot2::aes(fill = ggplot2::stat(count_transform(ggplot2::after_stat(count))))
                        )
                    } else if (geom[1] == "points") {
                        p <- p + scattermore::geom_scattermore(
                            alpha = 0.3, pointsize = 1, pixels = c(bins, bins)
                        )
                    } else if (geom[1] == "pointdensity") {
                        p <- p +
                            ggpointdensity::stat_pointdensity(
                                geom = GeomScattermore(),
                                ggplot2::aes(color = ggplot2::stat(count_transform(ggplot2::after_stat(density)))),
                                pointsize = 1.2,
                                pixels = c(bins, bins)
                            ) +
                            ggplot2::scale_color_continuous(type = "viridis")
                    }

                    p <- p +
                        ggpubr::theme_pubr() +
                        ggplot2::scale_fill_continuous(type = "viridis") +
                        ggplot2::theme(
                            legend.position = "none",
                            axis.title = ggplot2::element_text(size = 25)
                        )

                    # Add axis labels with cofactor if requested
                    if (axis_full_labels) {
                        p <- p +
                            ggplot2::xlab(paste0(marker_x, " [", transform_fun_name, "(z/", cofactor_x, ")]")) +
                            ggplot2::ylab(paste0(marker_y, " [", transform_fun_name, "(z/", cofactor_y, ")]"))
                    } else {
                        p <- p +
                            ggplot2::xlab(paste0(transform_fun_name, "(z/", cofactor_x, ")")) +
                            ggplot2::ylab(paste0(transform_fun_name, "(z/", cofactor_y, ")"))
                    }
                }

                if (length(add_ggplot_elements) != 0) {
                    for (element in add_ggplot_elements) {
                        p <- p + element
                    }
                }
                p_m <- list(p)
                names(p_m) <- xy_str

                plots_all_patchwork <- c(plots_all_patchwork, p_m)
                all_combos_str[[which_matching]] <- NULL
            } else {
                plots_all_patchwork <- c(plots_all_patchwork, list(patchwork::plot_spacer()))
            }
        }
    }

    # Construct full patchwork layout
    if (diag_plot) {
        plotlist_wrapper <- c(list(patchwork::plot_spacer()), xy_plots, plots_all_patchwork)
        ncol_wrapper <- length(marker_names) + 1
    } else {
        plots_all_patchwork[1:(length(marker_names) + 1)] <- c(list(patchwork::plot_spacer()), xy_plots)
        plotlist_wrapper <- plots_all_patchwork
        plotlist_wrapper[(seq_along(plotlist_wrapper)) %% (length(marker_names) + 1) == 0] <- NULL
        ncol_wrapper <- length(marker_names)
    }

    patchwork::wrap_plots(
        plotlist_wrapper,
        ncol = ncol_wrapper,
        byrow = TRUE,
        widths = c(0.15, rep(1, ncol_wrapper - 1)),
        heights = c(0.15, rep(1, ncol_wrapper - 1))
    )
}
