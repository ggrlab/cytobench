#' Plot Comparison of Feature Alignments Across Devices and Scaling Methods
#'
#' This function visualizes the density distributions of selected features across multiple devices and scaling methods.
#' It takes a nested list of data tables, applies specified scaling functions, computes density estimates, and plots
#' the results using ggplot2 with faceting by feature and scaling method.
#'
#' @param dt_listlist
#' A list of lists of data.tables. The outer list corresponds to the alignment method, the inner list device (which should be aligned). The innermost element is a data.table of feature measurements.
#' @param relevant_features Character vector of feature names to plot, referring to (some of) the column names of the data.tables. If NULL, all features present in the data are plotted.
#' @param n_resolution Integer specifying the number of points to use for density approximation (default: 1000).
#' @param scale_column_funs Named list of functions to scale feature values. Each function should accept a numeric vector and return a transformed numeric vector.
#'
#' @return A ggplot object showing density curves for each feature, device, and scaling method.
#'
#' @examples
#' # Example usage:
#' # plot_compare_alignment(dt_listlist, relevant_features = c("feature1", "feature2"))
#'
#' @export
plot_compare_alignment <- function(
    dt_listlist,
    relevant_features = NULL,
    n_resolution = 1e3,
    scale_column_funs = list("default" = function(x) {
        asinh(x * 100)
    }, "raw" = function(x) {
        asinh(x / 1e3)
    }), ylim = c(0, .4),
    faceting = c("horizontal", "vertical")) {
    rs_sn_long <- lapply(dt_listlist, data.table::rbindlist, idcol = "device", fill = TRUE)
    rescaled_variants <- data.table::rbindlist(rs_sn_long, idcol = "scale_column_fun", fill = TRUE)
    rescaled_variants_long <- data.table::melt(
        rescaled_variants,
        id.vars = c("device", "scale_column_fun"),
        variable.name = "feature",
        value.name = "value"
    )
    densfuns_ranges <- rescaled_variants_long[,
        {
            if (scale_column_fun %in% names(scale_column_funs)) {
                value <- scale_column_funs[[scale_column_fun]](value)
            } else {
                value <- scale_column_funs[["default"]](value)
            }
            resfun <- function(x) {
                NA
            }
            try(resfun <- approxfun(density(value, na.rm = TRUE)), silent = TRUE)
            .(
                dens_approxfun = list(list(resfun)),
                min = min(value),
                max = max(value)
            )
        },
        by = c("device", "scale_column_fun", "feature")
    ]
    approximated_dens_all <- data.table()
    for (feature_x in unique(densfuns_ranges$feature)) {
        current <- densfuns_ranges[feature == feature_x]
        range_total <- range(current$min, current$max, na.rm = TRUE)

        approximated_dens <- data.table()
        for (i in seq_along(current$dens_approxfun)) {
            afun <- current[["dens_approxfun"]][[i]][[1]]
            range_current <- c(range_total[1], range_total[2])
            x <- seq(range_current[[1]], range_current[[2]], length.out = n_resolution)
            approximated_dens <- rbind(
                approximated_dens, cbind(
                    current[i],
                    data.table(
                        x = x,
                        y = afun(x)
                    )
                )
            )
        }
        approximated_dens_all <- rbind(
            approximated_dens_all,
            approximated_dens
        )
    }
    approximated_dens_all[, scale_column_fun := factor(
        scale_column_fun,
        levels = unique(c("raw", "minmax", "relative", names(dt_listlist)))
    )]
    if (all(is.null(relevant_features))) {
        relevant_features <- unique(approximated_dens_all$feature)
    }
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
    if (!all(is.na(ylim))) {
        p0 <- p0 + ggplot2::ylim(ylim)
    }
    return(p0)
}
