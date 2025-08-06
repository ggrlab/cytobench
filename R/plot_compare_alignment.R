#' Plot Comparison of Feature Alignments Across Devices and Scaling Methods
#'
#' Visualizes the density distributions of selected features across multiple devices and scaling methods.
#' Each method-device-feature combination is plotted using an estimated density function. Useful to assess
#' alignment consistency or transformation effects across datasets.
#'
#' It takes a nested list of data tables, applies specified scaling functions, computes density estimates, and plots
#' the results using ggplot2 with faceting by feature and scaling method.
#' @param dt_listlist A named list of lists of `data.table`s. The outer list corresponds to the scaling or
#' alignment method (e.g. "raw", "minmax"). The inner list groups by device name (e.g., "navios", "cytoflex") - which should usually be aligned.
#' Each innermost element is a `data.table` of feature measurements.
#' @param relevant_features Character vector of feature names (column names) to include in the plot.
#' If `NULL`, all features found in the data are used.
#' @param n_resolution Integer. Number of x-values per density curve (default: 1000).
#' @param scale_column_funs Named list of transformation functions to apply to the feature values before
#' computing densities. Each function must take a numeric vector and return a transformed numeric vector.
#' @param ylim Numeric vector of length 2. Y-axis limits (e.g., c(0, 0.4)). Default: `c(0, 0.4)`.
#' @param faceting Either `"horizontal"` or `"vertical"`. Determines how ggplot facets are arranged.
#'
#' @return A `ggplot2` object showing overlaid density plots, faceted by feature and scaling method.
#'
#' @examples
#' # Assuming dt_listlist is a nested list of data.tables with aligned features:
#' # plot_compare_alignment(dt_listlist, relevant_features = c("FL1-A", "FL2-A"))
#'
#' @export
#' @keywords cytometry
#' @keywords relativisation
#' @examples
#' # Mock data: create 2 scaling methods, each with 2 devices, each with 3 features and 500 values
#' set.seed(123)
#' dt_listlist <- list(
#'     "raw" = list(
#'         "device1" = data.table(A = rnorm(500, 1), B = rnorm(500, 5), C = rnorm(500, 10)),
#'         "device2" = data.table(A = rnorm(500, 1.2), B = rnorm(500, 4.8), C = rnorm(500, 9.5))
#'     ),
#'     "default" = list(
#'         "device1" = data.table(A = rnorm(500, 1), B = rnorm(500, 5), C = rnorm(500, 10)),
#'         "device2" = data.table(A = rnorm(500, 1.1), B = rnorm(500, 5.1), C = rnorm(500, 10.5))
#'     )
#' )
#'
#' # Call the function with this mock data
#' p <- plot_compare_alignment(
#'     dt_listlist = dt_listlist,
#'     relevant_features = c("A", "B", "C"),
#'     ylim = NA,
#'     scale_column_funs = list("default" = function(x) {
#'         x
#'     })
#' )
#' print(p)
#'
#' fs <- simulate_fs(
#'     n_samples = 2,
#'     flowcore = FALSE,
#'     ncells = 250,
#'     columns = c("FL1-A", "FL2-A", "FL3-A")
#' )
#' fs_example <- list(
#'     "raw" = fs,
#'     "modified" = lapply(fs, function(x) {
#'         (x + 5) * 1.5
#'     })
#' )
#' p2 <- plot_compare_alignment(
#'     fs_example,
#'     ylim = NA,
#'     scale_column_funs = list(
#'         "default" = function(x) {
#'             x
#'         }
#'     )
#' )
#' print(p2)
#' p3 <- plot_compare_alignment(
#'     fs_example,
#'     ylim = NA,
#'     scale_column_funs = list(
#'         "raw" = function(x) {
#'             x
#'         },
#'         "modified" = function(x) {
#'             log10(x)
#'         }
#'     )
#' )
#' print(p3)
#'
plot_compare_alignment <- function(
    dt_listlist,
    relevant_features = NULL,
    n_resolution = 1e3,
    scale_column_funs = list(
        "default" = function(x) asinh(x * 100),
        "raw" = function(x) asinh(x / 1e3)
    ),
    ylim = c(0, .4),
    faceting = c("horizontal", "vertical")) {
    . <- feature <- x <- y <- device <- scale_column_fun <- NULL # to avoid R CMD check note about undefined global variable
    # Combine each method-device list into one data.table, add 'device' column
    rs_sn_long <- lapply(dt_listlist, data.table::rbindlist, idcol = "device", fill = TRUE)
    # Combine all scaling methods into one long data.table, add 'scale_column_fun'
    rescaled_variants <- data.table::rbindlist(rs_sn_long, idcol = "scale_column_fun", fill = TRUE)
    # Convert to long format: one row per (device, scaling method, feature value)
    rescaled_variants_long <- data.table::melt(
        rescaled_variants,
        id.vars = c("device", "scale_column_fun"),
        variable.name = "feature",
        value.name = "value"
    )

    # Compute one density approximator per (device, scaling method, feature)
    densfuns_ranges <- rescaled_variants_long[,
        {
            value <- if (scale_column_fun %in% names(scale_column_funs)) {
                scale_column_funs[[scale_column_fun]](value)
            } else {
                scale_column_funs[["default"]](value)
            }

            # Safely construct an approximated density function
            resfun <- function(x) NA
            try(resfun <- stats::approxfun(stats::density(value, na.rm = TRUE)), silent = TRUE)

            .(
                dens_approxfun = list(list(resfun)),
                min = min(value, na.rm = TRUE),
                max = max(value, na.rm = TRUE)
            )
        },
        by = c("device", "scale_column_fun", "feature")
    ]

    # Evaluate all density functions for each feature
    approximated_dens_all <- data.table()
    for (feature_x in unique(densfuns_ranges$feature)) {
        current <- densfuns_ranges[feature == feature_x]
        range_total <- range(current$min, current$max, na.rm = TRUE)

        for (i in seq_along(current$dens_approxfun)) {
            afun <- current$dens_approxfun[[i]][[1]]
            x_vals <- seq(range_total[1], range_total[2], length.out = n_resolution)

            dens_df <- data.table(
                x = x_vals,
                y = afun(x_vals)
            )
            # Add method and device labels
            dens_df <- cbind(current[i], dens_df)
            approximated_dens_all <- rbind(approximated_dens_all, dens_df)
        }
    }

    # Normalize factor levels for consistency in plots
    approximated_dens_all[, scale_column_fun := factor(
        scale_column_fun,
        levels = unique(c("raw", "minmax", "relative", names(dt_listlist)))
    )]

    # Auto-detect features if none provided
    if (all(is.null(relevant_features))) {
        relevant_features <- unique(approximated_dens_all$feature)
    }

    # Base ggplot object
    p0 <- ggplot2::ggplot(
        approximated_dens_all[feature %in% relevant_features],
        ggplot2::aes(x = x, y = y, col = device)
    ) +
        ggplot2::geom_line() +
        ggpubr::theme_pubr() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.line.y = ggplot2::element_blank()
        )

    # Add facet layout
    if (faceting[1] == "horizontal") {
        p0 <- p0 +
            ggh4x::facet_grid2(
                feature ~ scale_column_fun,
                scales = "free",
                independent = "all"
            )
    } else {
        p0 <- p0 +
            ggh4x::facet_grid2(
                scale_column_fun ~ feature,
                scales = "free",
                independent = "all"
            )
    }

    # Optional y-limit
    if (!all(is.na(ylim))) {
        p0 <- p0 + ggplot2::ylim(ylim)
    }

    return(p0)
}
