extract_quantile <- function(fcs_dir = "data-raw/s001",
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
            return(NULL)
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
            stop("Not implemented")
            extract_relevant_mfis_multistain(
                loaded_fcs_multistain,
                transform_fun = transform_fun,
                relevant_columns = multistain_columns
            )
        },
        error = function(e) {
            return(NULL)
        }
    )
    # Unname each column
    return(
        list(
            "single" = joint_df,
            "multi" = relevant_mfis_multi
        )
    )
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
    approxmethod = function(...) {
        # rule = 2 ensures that the function can extrapolate outside the range of the data
        approxfun(..., rule = 2)
    },
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
            ff_exprs <- flowCore::exprs(ff)
            expr_vals <- ff_exprs[, target_channel]
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
                    if (!marker %in% colnames(ff_exprs)) {
                        return(NULL)
                    }
                    y_vals <- ff_exprs[cluster_labels == "positive", marker]
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
