#' @title Adaptive Break Generator for ggplot2 Scales
#'
#' @description
#' Generates aesthetically pleasing axis breaks that automatically switch
#' between *linear (pretty)* and *logarithmic-like* spacing depending on the
#' data range. This is particularly useful for plotting transformed data
#' (e.g., `asinh`-scaled cytometry values) where ranges can vary widely.
#'
#' @param x `numeric()`
#'   Numeric vector from which to determine the axis range and break placement.
#' @param n `integer(1)`
#'   Desired approximate number of breaks to generate. Defaults to `6`.
#'
#' @details
#' The function inspects the numeric range of `x`:
#' - If the total span (`max(x) - min(x)`) is below 10, it uses
#'   [scales::breaks_pretty()] to produce human-readable, evenly spaced ticks.
#' - Otherwise, it applies [scales::breaks_log()] to the absolute values of `x`
#'   and restores the original sign to obtain symmetric log-like breaks that
#'   handle both positive and negative ranges gracefully.
#'
#' This behavior allows for consistent axis labeling across datasets with
#' highly variable magnitude, while maintaining interpretability near zero.
#'
#' @return
#' A numeric vector of break positions suitable for use in ggplot2 scale
#' functions (e.g., `scale_x_continuous(breaks = breaks_auto)`).
#'
#' @examples
#' x <- c(-1000, -10, 0, 10, 1000)
#' breaks_auto(x)
#'
#' ggplot(data.frame(x = -100:100, y = -100:100), aes(x, y)) +
#'     geom_point() +
#'     scale_x_continuous(breaks = breaks_auto) +
#'     scale_y_continuous(breaks = breaks_auto)
#'
#' @export
breaks_auto <- function(x, n = 5) {
    # Compute data range, ignoring missing values
    rng <- range(x, na.rm = TRUE)

    # If the range is small, use evenly spaced "pretty" breaks
    if (diff(rng) < 1e4) {
        scales::breaks_pretty(n = n)(x)
    } else {
        scales::breaks_log(n = n)(abs(x))
    }
}


#' Plot simple pre-/post-gate 2D scatter views from a GatingHierarchy
#'
#' @description
#' Given a single-sample \strong{GatingSet} or a \strong{GatingHierarchy}, this helper
#' draws a 2D scatter for each population, showing the parent (pre-gate) and the
#' population data (post-gate). Axes are transformed according to the gate's
#' channel transforms obtained from `flowWorkspace::gh_get_transformations()`.
#'
#' TODO:
#'  It would be great to have the same parent gate (and x-y combination (?) plot the cells from
#'  all its children in different colors. E.g. for a a bivariate plot with 3 rectangle gates,
#'  plot the pregated only once, and the 3 gated populations in different colors.
#' @param gatingset A single-sample `flowWorkspace::GatingSet` \emph{or}
#'   a `flowWorkspace::GatingHierarchy`. If a `GatingSet` is provided, it must
#'   contain exactly one sample.
#' @param populations `NULL` or a character vector of population paths to plot.
#'   If `NULL` (default), all non-root populations are used.
#' @param gg_layer A ggplot2 layer (e.g. `scattermore::geom_scattermore(alpha = 0.4)`
#'   or `ggplot2::geom_point(alpha = 0.2, size = 0.2)`) used to draw points.
#'   Defaults to a fast scatter layer from **scattermore**.
#' @param facet Logical, if `TRUE` (default) facets columns by `gatingstatus`
#'   ("pregate" vs "postgate"). If `FALSE`, returns a single panel with both
#'   statuses overplotted
#' @param verbose Logical, suppresses informative messages for non-2D gates or
#'   missing transforms. Default `FALSE`.
#' @param ... Unused; reserved for future extensions.
#'
#' @return A named `list` of `ggplot` objects, named by population path.
#'
#' @details
#' \strong{On transforms:} ggplot2 expects a \emph{trans object} from **scales**;
#' the correct constructor is `scales::trans_new(name, transform, inverse, domain)`.
#' This function builds such trans objects from the forward and inverse functions
#' returned by `gh_get_transformations()`. Populations whose channels lack either
#' forward or inverse transforms will be plotted on the untransformed scale with
#' a warning (unless `verbose = TRUE`).
#'
#' \strong{Limitations:} Only 2-D gates are supported (typical rectangle/polygon
#' gates). Non-2-D populations are skipped with an informative message.
#'
#' @examples
#' \dontrun{
#' # Using scattermore for fast plotting
#' plots <- plot_gating_simple(
#'     gatingset = my_single_sample_gs,
#'     populations = c("/Cells", "/Cells/Singlets", "/Cells/Singlets/Live")
#' )
#' # Print one:
#' print(plots[["/Cells/Singlets"]])
#'
#' # Use base geom_point and overplot pre/post in one panel
#' plots2 <- plot_gating_simple(
#'     my_gh,
#'     gg_layer = ggplot2::geom_point(alpha = 0.15, size = 0.2),
#'     facet = FALSE
#' )
#' }
#'
#' @export
plot_gating_simple <- function(
    gatingset,
    transformlist = NULL,
    populations = NULL,
    gg_layer = scattermore::geom_scattermore(alpha = 0.4),
    facet = TRUE,
    verbose = FALSE,
    ...) {
    # ---- Normalize input to a single GatingHierarchy ----
    if (inherits(gatingset, "GatingHierarchy")) {
        gh <- gatingset
    } else if (inherits(gatingset, "GatingSet")) {
        # Enforce single-sample GS
        if (length(gatingset) != 1L) {
            stop("Please provide a single-sample GatingSet, or a single GatingHierarchy.")
        }
        gh <- gatingset[[1]]
    } else {
        stop("`gatingset` must be a GatingHierarchy or a single-sample GatingSet.")
    }

    # ---- Obtain forward and inverse transformations per channel ----
    gh_copy <- flowWorkspace::gs_clone(gh)[[1]]
    gh_copy <- flowWorkspace::transform(gh_copy, transformlist)
    transforms_fwd <- flowWorkspace::gh_get_transformations(gh_copy)
    transforms_inv <- flowWorkspace::gh_get_transformations(gh_copy, inverse = TRUE)
    rm(gh_copy)
    trans_map <- sapply(
        names(transforms_fwd),
        simplify = FALSE,
        function(fun_x) {
            # # Create a new scale function that applies the transformation
            scales::trans_new(
                name = fun_x,
                transform = transforms_fwd[[fun_x]],
                inverse = transforms_inv[[fun_x]],
                breaks = breaks_auto
            )
        }
    )

    # ---- Determine which populations to plot ----
    if (!is.null(populations)) {
        paths <- populations
    } else {
        paths <- flowWorkspace::gh_get_pop_paths(gh)[-1] # drop root
    }

    # ---- Helper: build one plot for a given population path ----
    build_plot <- function(path_x) {
        # Get gate and parent; skip if unavailable
        current_gate <- try(flowWorkspace::gh_pop_get_gate(gh, path_x), silent = TRUE)
        if (inherits(current_gate, "try-error")) {
            if (!verbose) message("Skipping '", path_x, "': cannot retrieve gate.")
            return(NULL)
        }

        parent_path <- try(flowWorkspace::gh_pop_get_parent(gh, path_x), silent = TRUE)
        if (inherits(parent_path, "try-error")) {
            if (!verbose) message("Skipping '", path_x, "': cannot retrieve parent population.")
            return(NULL)
        }

        gate_params <- current_gate@parameters
        if (length(gate_params) != 2L) {
            if (!verbose) {
                message(path_x, "': gate is not 2D (found ", length(gate_params), " params). Using the first column as second axis.")
            }
            gate_params <- c(gate_params, flowWorkspace::colnames(gh)[1])
        }
        x_y <- names(gate_params)

        # Realize pre- and post-gate data views and extract the two columns
        dt_pregate <- flowWorkspace::gh_pop_get_data(gh, parent_path)[, x_y, drop = FALSE] |>
            flowWorkspace::realize_view() |>
            flowCore::exprs() |>
            data.table::as.data.table()
        dt_postgate <- flowWorkspace::gh_pop_get_data(gh, path_x)[, x_y, drop = FALSE] |>
            flowWorkspace::realize_view() |>
            flowCore::exprs() |>
            data.table::as.data.table()

        # Bind rows with status indicator
        dt_joint <- data.table::rbindlist(
            list("pregate" = dt_pregate, "postgate" = dt_postgate),
            idcol = "gatingstatus"
        ) |>
            dplyr::mutate(gatingstatus = factor(gatingstatus, levels = c("pregate", "postgate")))

        # Prepare axis transforms if available; warn if missing
        x_trans <- trans_map[[x_y[1]]]
        y_trans <- trans_map[[x_y[2]]]
        if (is.null(x_trans) && !verbose) {
            message("Path '", path_x, "': no transform (or inverse) for channel '", x_y[1], "'. Using identity.")
            x_trans <- scales::identity_trans()
        }
        if (is.null(y_trans) && !verbose) {
            message("Path '", path_x, "': no transform (or inverse) for channel '", x_y[2], "'. Using identity.")
            y_trans <- scales::identity_trans()
        }

        # Build plot
        p <- ggplot2::ggplot(
            dt_joint,
            ggplot2::aes(x = .data[[x_y[1]]], y = .data[[x_y[2]]])
        ) +
            gg_layer + # user-provided/ default point layer
            ggpubr::theme_pubr() + # clean theme
            ggplot2::scale_x_continuous(trans = x_trans) +
            ggplot2::scale_y_continuous(trans = y_trans) +
            ggplot2::ggtitle(path_x)

        # Facet or overplot depending on user choice
        if (isTRUE(facet)) {
            p <- p + ggplot2::facet_grid(. ~ gatingstatus)
        }

        p
    }

    # ---- Build plots for all requested paths ----
    out <- sapply(paths, simplify = FALSE, build_plot)

    # Drop NULLs (unplottable nodes) but keep names on remaining
    keep <- !vapply(out, is.null, logical(1))
    if (!any(keep)) {
        if (!verbose) message("No plottable 2D populations were found.")
        return(list())
    }

    out[keep]
}
