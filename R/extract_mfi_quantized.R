extract_mfi_quantized <- function(fcs_dir = "data-raw/s001",
                                  regex_singlestain = "(-(CD3-.*)|(none))\\.fcs$",
                                  transform_fun = function(x) {
                                      asinh(x / 1e3)
                                  },
                                  regex_multistain = "(_15-MasterMix)\\.fcs",
                                  multistain_columns = c(
                                      "FITC-A", "PE-A", "ECD-A", "PC5.5-A", "PC7-A",
                                      "APC-A", "AF700-A", "AA750-A", "PB-A", "KrO-A"
                                  ),
                                  gating_set_file = NULL,
                                  gate_extract = NULL,
                                  ...) {
    browser()
    joint_df <- tryCatch(
        {
            loaded_fcs <- cytobench:::load_mfi_files(
                fcs_dir = fcs_dir,
                regex_singlestain = regex_singlestain,
                gating_set_file = gating_set_file,
                gate_extract = gate_extract
            )
            # If multistaining is enabled, the following extracts the potential UNSTAINED sample
            extract_singlestain_quantilefun(loaded_fcs, transform_fun = transform_fun, relevant_columns = multistain_columns)
        },
        error = function(e) {
            tibble::tibble(
                "feature" = NA,
                "negative" = NA,
                "positive" = NA,
                "unstained" = NA,
                "negative.sd" = NA,
                "positive.sd" = NA,
                "unstained.sd" = NA
            )
        }
    )

    relevant_mfis_multi <- tryCatch(
        {
            # Then extract the actually multi-stained sample(s)
            loaded_fcs_multistain <- load_mfi_files(
                fcs_dir = fcs_dir,
                regex_singlestain = regex_multistain,
                gating_set_file = gating_set_file,
                gate_extract = gate_extract
            )
            extract_relevant_mfis_multistain(
                loaded_fcs_multistain,
                transform_fun = transform_fun,
                relevant_columns = multistain_columns
            )
        },
        error = function(e) {
            return(tibble::tibble(
                "sample" = NA,
                "feature" = NA,
                "negative" = NA,
                "positive" = NA,
                "unstained" = NA,
                "negative.sd" = NA,
                "positive.sd" = NA,
                "unstained.sd" = NA
            ))
        }
    )
    if (nrow(relevant_mfis_multi) > 0) {
        #  Merge single and multi stainings
        joint_df <- dplyr::left_join(
            joint_df,
            relevant_mfis_multi,
            by = "feature",
            suffix = c("", ".multi")
        )
    }
    # The following replaces the values of the respective columns
    # by their multi-stained counterparts if they are NA.
    for (x in c("positive", "negative", "positive.sd", "negative.sd")) {
        if (all(is.na(joint_df[[x]]))) {
            joint_df[[x]] <- joint_df[[paste0(x, ".multi")]]
        }
    }

    # Unname each column
    return(dplyr::mutate(joint_df, dplyr::across(tidyr::everything(), ~ unname(.))))
}

extract_singlestain_quantilefun <- function(
    loaded_fcs,
    n_quantiles = 100,
    quantiles = NULL,
    transform_fun = function(x) {
        asinh(x / 1e3)
    },
    approxmethod = approxfun,
    relevant_columns,
    ...) {
    if (is.null(quantiles)) {
        # Create default quantiles if not provided
        quantiles <- seq(from = 0, to = 1, length.out = n_quantiles + 2)
        # Remove the first and last quantiles to avoid exact 0 and 1
        quantiles <- quantiles[-c(1, length(quantiles))]
    }


    qfuns <- sapply(names(loaded_fcs), simplify = FALSE, function(f_x) {
        ff_x <- flowWorkspace::cytoframe_to_flowFrame(loaded_fcs[[f_x]][[1]])
        nonempty_channel <- which(flowCore::markernames(ff_x) != "empty")
        nonempty_name <- names(nonempty_channel)

        # Check if there is more than one non-empty channel
        if (length(nonempty_channel) > 1) {
            warning("More than one non-empty channel in ", f_x)
            stop("Untested")
        }

        # If no non-empty channel, it is the unstained sample
        if (length(nonempty_channel) == 0) {
            # stop("Untested")
            return(NULL)
        } else {
            # Cluster the relevant channel into two populations
            values_nonempty <- flowCore::exprs(ff_x)[, nonempty_name]
            clustering <- stats::kmeans(transform_fun(values_nonempty), centers = 2)
            positive_cluster <- which.max(clustering$centers)
            clusters_named <- ifelse(clustering$cluster == positive_cluster, "positive", "negative")
            # per cluster, create an approximate function where x=quantile and y=MFI
            ff_dt <- data.table::as.data.table(flowCore::exprs(ff_x))
            ff_dt[, cluster := clusters_named]
            ff_dt_positive <- ff_dt[cluster == "positive"]
            # The following is just for testing
            # tmp <- stats::ecdf(ff_dt_positive[["FITC-A"]])
            # ff_dt_positive[, ecdf := tmp(`FITC-A`)]
            # ff_dt_positive[, ecdf_quantiled := quantile(ff_dt_positive[["FITC-A"]], ecdf)]
            # summary(ff_dt_positive[["ecdf_quantiled"]] - ff_dt_positive[["FITC-A"]])
            # summary(ff_dt_positive[["ecdf_quantiled"]] / ff_dt_positive[["FITC-A"]])
            # plot(ff_dt_positive[["ecdf_quantiled"]], ff_dt_positive[["FITC-A"]], type = "l")
            # abline(0, 1, col = "red")
            # data.table::setkeyv(ff_dt_positive, "ecdf")

            # pdf("removeme2.pdf")
            # ggplot(
            #     ff_dt_positive |>
            #         tidyr::pivot_longer(
            #             cols = -c("ecdf", "cluster"),
            #             names_to = "marker",
            #             values_to = "value"
            #         ) |> dplyr::filter(!marker %in% c("TIME", "SSC-A", "FSC-A", "FSC-W", "FITC-A")), aes(x = ecdf, y = value, col = marker)
            # ) +
            #     geom_line() +
            #     ggpubr::theme_pubr() +
            #     facet_wrap(~marker, scales = "free") +
            #     ggplot2::geom_smooth(formula = y ~ x, method = "loess", col = "black")
            # dev.off()
            approxfun_per_marker <- apply(ff_dt_positive[, -ncol(ff_dt_positive), with = FALSE], 2, function(x) {
                # Note: Do NOT use ff_dt completely here to also get the negative populations.
                # You should _probably_ calculate the values corresponding to the qunatiles of "THE" positive population,
                # not based on the NEW quantiles from each negative population

                # Create an approximation function for the quantile-MFI relationship
                # This uses the specified approximation method (e.g., splinefun or approxfun)
                # x: quantiles, y: quantile values of the channel, ...: additional arguments
                approximation_fun <- approxmethod(
                    x = quantiles,
                    y = quantile(x, quantiles, names = FALSE),
                    ...
                )
                return(approximation_fun)
            })
            approxfun_nonempty <- approxfun_per_marker[nonempty_name]
        }

        return(approxfun_nonempty)
    })

    # The following shows the generated quantile functions (x=quantile, y=respective MFI)
    # and we observe that the "Quantile-MFI" actually deviates considerably from the
    # previously used MFI.
    # Our base relativisation approach would have "taken" the median fluorescence intensity
    # of the whole positive population, (the 50% quantile).
    # plotting_df <- lapply(names(qfuns_unlisted), function(x) {
    #     data.table::data.table(
    #         marker = x,
    #         quantiles = quantiles,
    #         approximation = qfuns_unlisted[[x]](quantiles)
    #     ) |> tidyr::pivot_longer(
    #         cols = -c("quantiles", "marker"),
    #         names_to = "cluster",
    #         values_to = "value"
    #     )
    # }) |>
    #     data.table::rbindlist() |>
    #     dplyr::mutate(
    #         pos_neg = stringr::str_extract(marker, "positive|negative"),
    #         markername = stringr::str_extract(marker, "\\..*$")
    #     )
    # pdf("removeme.pdf", width = 12)
    # ggplot(plotting_df |> dplyr::filter(
    #     !grepl("TIME", marker),
    #     !grepl("SC", marker),
    # ), aes(x = quantiles, y = value, col = markername, group = marker)) +
    #     geom_line() +
    #     scale_y_log10()
    # dev.off()
    # Extract median fluorescence intensities (MFIs) from the loaded FCS files
    return(qfuns)
}


#' Extract Quantile Approximation Functions for Single-Stain Flow Cytometry Data
#'
#' This function processes a list of loaded FCS files containing single-stain samples,
#' identifies the relevant marker channel, applies a transformation to the expression values,
#' separates positive and negative populations using k-means clustering, and computes
#' quantile approximation functions for the positive population.
#'
#' @param loaded_fcs
#' A named list of loaded FCS objects, where each element is a cytoset. Usually created by
#' \code{cytobench:::load_mfi_files}. 
#' @param n_quantiles Integer.
#' Number of quantiles to use if \code{quantiles} is not provided. Default is 100.
#' @param quantiles Numeric vector.
#'  Specific quantile probabilities to use (between 0 and 1). If \code{NULL}, defaults to evenly spaced quantiles.
#' @param transform_fun Function.
#' Transformation to apply to expression values. Only used for clustering the (then) transformed values.
#' Default is \code{asinh(x / 1e3)}.
#' @param approxmethod Function.
#' Method to approximate the quantile function (e.g., \code{approxfun}). Default is \code{approxfun}.
#' @param relevant_columns Character vector.
#' Names of marker columns to compute quantiles for.
#' @param ...
#' Additional arguments passed to the approximation method.
#'
#' @return
#' A named list of lists, where each element corresponds to a sample and contains a quantile approximation
#'  function for the relevant marker channel.
#'
#' @details
#' For each sample, the function:
#' \itemize{
#'   \item Identifies the non-empty marker channel.
#'   \item Applies the specified transformation to the expression values.
#'   \item Uses k-means clustering to separate positive and negative populations.
#'   \item Restricts analysis to the positive population.
#'   \item Computes quantile approximation functions for the specified marker columns.
#' }
#' Samples with no non-empty channels are skipped. Samples with multiple non-empty channels result in an error.
#'
#'
#' @export
extract_singlestain_quantilefun <- function(
    loaded_fcs,
    n_quantiles = 100,
    quantiles = NULL,
    transform_fun = function(x) asinh(x / 1e3),
    approxmethod = approxfun,
    relevant_columns,
    ...) {
    # Create default quantiles if not provided
    if (is.null(quantiles)) {
        quantiles <- seq(0, 1, length.out = n_quantiles + 2)[-c(1, n_quantiles + 2)]
    }

    qfuns <- sapply(
        names(loaded_fcs),
        simplify = FALSE,
        function(f_x) {
            cyto_obj <- loaded_fcs[[f_x]][[1]] # I expect a cytoset with each 1 sample.
            ff <- flowWorkspace::cytoframe_to_flowFrame(cyto_obj)
            marker_names <- flowCore::markernames(ff)
            nonempty_idx <- which(marker_names != "empty")
            nonempty_names <- names(nonempty_idx)

            # If no non-empty channel: skip (unstained sample)
            if (length(nonempty_idx) == 0) {
                return(NULL)
            }

            # If multiple non-empty channels: error for now
            if (length(nonempty_idx) > 1) {
                stop(sprintf("Sample %s has multiple non-empty channels. Currently unsupported.", f_x))
            }

            target_channel <- nonempty_names

            # Extract and transform expression values
            expr_vals <- flowCore::exprs(ff)[, target_channel]
            transformed_vals <- transform_fun(expr_vals)

            # Perform k-means clustering to separate positive and negative populations
            clustering <- stats::kmeans(transformed_vals, centers = 2)
            positive_cluster <- which.max(clustering$centers)
            cluster_labels <- ifelse(clustering$cluster == positive_cluster, "positive", "negative")

            # Compute per-marker quantile approximations
            qfun_list <- sapply(
                relevant_columns,
                simplify = FALSE,
                function(marker) {
                    if (!marker %in% colnames(ff)) {
                        return(NULL)
                    }
                    y_vals <- flowCore::exprs(ff)[cluster_labels == "positive", marker]
                    # Note: Do NOT use ff completely here to also get the negative populations.
                    # You should _probably_ calculate the values corresponding to the quantiles of "THE" positive population,
                    # not based on the NEW quantiles from each negative population

                    # Create an approximation function for the quantile-MFI relationship
                    # This uses the specified approximation method (e.g., splinefun or approxfun)
                    # x: quantiles, y: quantile values of the channel, ...: additional arguments
                    approxmethod(
                        x = quantiles,
                        y = stats::quantile(y_vals, probs = quantiles, names = FALSE),
                        ...
                    )
                }
            )
            return(qfun_list[target_channel])
        }
    )
    return(qfuns)
}
