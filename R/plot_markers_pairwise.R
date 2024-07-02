#' Plot all pairwise marker combinations with given cofactors
#'
#' Plots all pairwise marker combinations with given cofactors. The cofactors are used to scale the markers before plotting. The plot is a hexbin plot of the asinh transformed values of the markers. The cofactors are used to scale the markers before plotting.
#' @param ff A flowFrame object
#' @param cofactor_namedvec A named vector with the cofactors for each marker. The names of the vector should be the same as the marker names in the flowFrame object.
#' @param special_cofactor_list A named list with special cofactors for certain marker combinations. The names of the list should be the marker combinations, separated by an underscore. The values should be a vector of two values, the first for the x-axis marker and the second for the y-axis marker.
#' @param bins The number of bins for the hexbin plot
#' @param diag_plot Whether to include the diagonal plots
#' @param debugplots Whether to include debug plots (No cells, only text to show the layout)
#' @param axis_full_labels Whether to include the cofactors in the axis labels
#' @param n_cells The number of cells to plot. Default is Inf, which plots all cells.
#' @export
#' @examples
#' pdf("removeme.pdf", width = 60, height = 50)
#' print(
#'     plot_combinations(
#'         current_sample_gated[[1]],
#'         pre_rescale_asinh_sscheck[[device_x]],
#'         special_asinh[[device_x]],
#'         debugplots = TRUE
#'     )
#' )
#' print(
#'     plot_combinations(
#'         current_sample_gated[[1]],
#'         pre_rescale_asinh_sscheck[[device_x]],
#'         special_asinh[[device_x]],
#'         diag_plot = FALSE,
#'         debugplots = TRUE
#'     )
#' )
#' dev.off()
#'
plot_markers_pairwise <- function(ff,
                              cofactor_namedvec,
                              special_cofactor_list,
                              bins = 50,
                              diag_plot = FALSE,
                              debugplots = FALSE,
                              axis_full_labels = TRUE,
                              n_cells = Inf) {
    xy_plots <- list()
    xy_plots_rotated <- list()
    for (marker_x in names(flowCore::markernames(ff))) {
        xy_plots[[marker_x]] <- ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0, y = 0, size = 18, label = marker_x) +
            ggplot2::theme_void()
        xy_plots_rotated[[marker_x]] <- ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0, y = 0, size = 18, label = marker_x, angle = 90) +
            ggplot2::theme_void()
    }
    # pdf("removeme.pdf")
    # xy_plots
    # dev.off()
    # # flowCore::markernames(gated)
    all_marker_combinations <- combn(
        names(flowCore::markernames(ff)),
        2,
        simplify = FALSE
    )

    all_marker_combinations_str <- lapply(all_marker_combinations, function(x) paste0(x, collapse = "_"))
    # plot_all_combinations <- list()
    plots_all_patchwork <- list()
    gated_exprs <- data.table::as.data.table(flowCore::exprs(ff))
    if (!is.infinite(n_cells)) {
        # https://stackoverflow.com/questions/24685421/how-do-you-extract-a-few-random-rows-from-a-data-table-on-the-fly
        gated_exprs <- gated_exprs[sample(.N, n_cells)]
    }
    for (marker_y in names(xy_plots)) {
        plots_all_patchwork <- c(plots_all_patchwork, xy_plots_rotated[marker_y])
        for (marker_x in names(xy_plots)) {
            xy_str <- paste0(marker_x, "_", marker_y)
            which_matching <- which(unlist(all_marker_combinations_str) == xy_str)
            # if(xy_str == "PE-A_ECD-A"){stop()}
            # if(xy_str == "PE-A_FITC-A"){stop()}
            if (xy_str %in% all_marker_combinations_str) {
                cofactor_x <- cofactor_namedvec[[marker_x]]
                cofactor_y <- cofactor_namedvec[[marker_y]]
                if (xy_str %in% names(special_cofactor_list)) {
                    cofactor_x <- special_cofactor_list[[xy_str]][1]
                    cofactor_y <- special_cofactor_list[[xy_str]][2]
                }

                if (debugplots) {
                    p_markers <- ggplot2::ggplot() +
                        ggplot2::annotate("text", x = 0, y = 0, size = 15, label = paste0(xy_str, "\nx: ", cofactor_x, "   : ", cofactor_y)) +
                        ggplot2::theme_void()
                } else {
                    dt_transformed <- tibble::tibble(
                        x = asinh(gated_exprs[[marker_x]] / cofactor_x),
                        y = asinh(gated_exprs[[marker_y]] / cofactor_y)
                    )
                    p_markers <- ggplot2::ggplot(dt_transformed, ggplot2::aes(x = x, y = y)) +
                        ggplot2::geom_hex(bins = bins) +
                        # geom_density_2d_filled() +
                        ggpubr::theme_pubr() +
                        ggplot2::scale_fill_continuous(type = "viridis", trans = "log10") +
                        ggplot2::theme(
                            legend.position = "none",
                            axis.title = ggplot2::element_text(size = 25)
                        )
                    if (axis_full_labels) {
                        p_markers <- p_markers +
                            ggplot2::xlab(paste0(marker_x, " [asinh(z/", cofactor_x, ")]")) +
                            ggplot2::ylab(paste0(marker_y, " [asinh(z/", cofactor_y, ")]"))
                    } else {
                        p_markers <- p_markers +
                            ggplot2::xlab(paste0("asinh(z/", cofactor_x, ")")) +
                            ggplot2::ylab(paste0("asinh(z/", cofactor_y, ")"))
                    }
                }
                p_m <- list(p_markers)
                names(p_m) <- xy_str
                plots_all_patchwork <- c(plots_all_patchwork, p_m)
                all_marker_combinations_str[[which_matching]] <- NULL
            } else {
                plots_all_patchwork <- c(plots_all_patchwork, list(patchwork::plot_spacer()))
            }
        }
    }
    if (diag_plot) {
        plotlist_wrapper <- c(list(patchwork::plot_spacer()), xy_plots, plots_all_patchwork)
        ncol_wrapper <- length(xy_plots) + 1
    } else {
        plots_all_patchwork[1:(length(xy_plots) + 1)] <- c(list(patchwork::plot_spacer()), xy_plots)
        plotlist_wrapper <- plots_all_patchwork
        plotlist_wrapper[(1:length(plotlist_wrapper)) %% (length(xy_plots) + 1) == 0] <- NULL
        ncol_wrapper <- length(xy_plots)
    }
    plots_wrapped <- patchwork::wrap_plots(
        plotlist_wrapper,
        ncol = ncol_wrapper,
        byrow = TRUE,
        widths = c(.15, rep(1, ncol_wrapper - 1)),
        heights = c(.15, rep(1, ncol_wrapper - 1))
    )
    return(plots_wrapped)
}
