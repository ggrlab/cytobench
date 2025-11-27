#' Plot all pairwise marker combinations with given cofactors
#'
#' Generates a matrix of 2D-density plots for all marker combinations in a flow cytometry dataset.
#' Axes can be transformed and scaled using marker-specific cofactors. This is particularly useful
#' for visual QC of compensation, cofactor choices and signal separation.
#'
#' @param ff_df A `flowFrame` or data.frame like  element containing the cytometry data.
#' @param cofactor_namedvec
#' A named numeric vector providing default cofactors per marker. If not given, defaults to 1500 for all column names in `ff_df`.
#' @param special_cofactor_list
#' A named list with overrides for specific marker combinations.
#' Each name should be `"markerX_markerY"`, with values being a
#' numeric vector of length 2: `c(cofactor_x, cofactor_y)`.
#' @param transform_fun A function used to transform intensity values (default: `asinh`).
#' @param transform_fun_name Character label for the transformation function (used in axis labels).
#' @param verbose Logical. If TRUE, prints progress messages to the console.
#' @param engine Character. One of `"ggplot"` (default) or `"base"` - determines the plotting engine used.
#' @param geom Character.
#' One of `"hex"` (default), `"points"`, or `"pointdensity"` -
#' determines the geometry used.
#' Speed concerns:
#'  - For engine ggplot2, hex is usually the fastest.
#'  - Fastest overall is `"points"` with engine `"base"`.
#' @param n_cells
#' Number of cells/rows to downsample randomly (default: `Inf` = use all).
#' @param count_transform
#' Function to transform hexbin/pointdensity counts or density estimates (e.g. `log10(count + 1)`).
#' @param title_global
#' Character. An optional global title for the entire plot.
#' @param modelines
#' Logical. If TRUE, adds red lines at the modes of each marker in the scatterplots.
#' @param bins `engine = "ggplot"`. Number of bins (resolution) in each 2D plot (default: 50). (hex geom)
#' @param diag_plot
#' Logical, `engine = "ggplot"`. Whether to include labeled marker names along the diagonal (default: FALSE).
#' @param debugplots
#' Logical, `engine = "ggplot"`. If TRUE, plots text placeholders instead of actual data - useful for layout testing.
#' @param axis_full_labels
#' Logical, `engine = "ggplot"`. If TRUE, axis labels include marker names and cofactor formulas.
#' @param add_ggplot_elements
#' `engine = "ggplot"`. A list of ggplot layers to add to each panel (e.g. `list(ggplot2::theme_minimal())`).
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
plot_markers_pairwise <- function(
    # Arguments for both cases
    ff_df,
    cofactor_namedvec,
    special_cofactor_list = list(),
    transform_fun = asinh,
    transform_fun_name = "asinh",
    verbose = FALSE,
    engine = c("ggplot", "base"),
    n_cells = Inf,
    title_global = NULL,
    modelines = TRUE,
    # Arguments specific to ggplot engine
    geom = c("hex", "points", "pointdensity"),
    bins = 50,
    diag_plot = FALSE,
    debugplots = FALSE,
    axis_full_labels = TRUE,
    count_transform = function(x) log10(x + 1),
    add_ggplot_elements = list(),
    ...) {
    if (missing(cofactor_namedvec)) {
        cofactor_namedvec <- setNames(
            rep(1500, ncol(ff_df)),
            flowCore::colnames(ff_df) # works for data.frames _and_ flowFrames
        )
    }
    if ("flowFrame" %in% class(ff_df)) {
        exprs_dt <- data.table::as.data.table(flowCore::exprs(ff_df)[, names(cofactor_namedvec)])
    } else {
        exprs_dt <- ff_df
        if (data.table::is.data.table(exprs_dt)) {
            exprs_dt <- exprs_dt[, names(cofactor_namedvec), with = FALSE]
        } else {
            exprs_dt <- data.table::data.table(exprs_dt[, names(cofactor_namedvec)])
        }
    }


    if (!is.infinite(n_cells)) {
        exprs_dt <- exprs_dt[sample(.N, min(.N, n_cells))]
    }

    if (engine[1] == "base") {
        plot_markers_pairwise_base(
            df = exprs_dt,
            cofactor_namedvec = cofactor_namedvec,
            special_cofactor_list = special_cofactor_list,
            transform_fun = transform_fun,
            transform_fun_name = transform_fun_name,
            geom = geom,
            title_global = title_global,
            modelines = modelines,
            ...
        )
    } else if (engine[1] == "ggplot") {
        plot_markers_pairwise_ggplot(
            df = exprs_dt,
            cofactor_namedvec = cofactor_namedvec,
            special_cofactor_list = special_cofactor_list,
            transform_fun = transform_fun,
            transform_fun_name = transform_fun_name,
            geom = geom,
            bins = bins,
            diag_plot = diag_plot,
            debugplots = debugplots,
            axis_full_labels = axis_full_labels,
            count_transform = count_transform,
            add_ggplot_elements = add_ggplot_elements,
            title_global = title_global,
            modelines = modelines,
            ...
        )
    }
}
