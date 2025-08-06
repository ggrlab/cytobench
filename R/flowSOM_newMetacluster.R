#' Predict and Reassign MetaClusters Based on a FlowSOM Result
#'
#' This function reassigns or updates meta-cluster labels for either:
#' (1) individual cells with known cluster assignments, or
#' (2) pre-aggregated cluster-level counts across samples.
#'
#' Optionally, a new number of meta-clusters can be specified, triggering recomputation via
#' `ConsensusClusterPlus`. The function also optionally updates the `FlowSOM` object to reflect
#' the new meta-clustering.
#'
#' @param flowsom_result A `FlowSOM` result object, as returned by `flowSOM()` or `flowSOM_optimal()`.
#' @param clustered_df A `data.frame` that either:
#'   \itemize{
#'     \item Contains a `"cluster"` column mapping each cell to a cluster (e.g., from a previous FlowSOM run), or
#'     \item Contains wide-format columns `"cluster_1"`, `"cluster_2"`, ..., representing counts of cells per cluster.
#'   }
#'   If provided, the function returns an updated version with newly assigned meta-clusters.
#' @param n_metacluster Optional. If provided, the FlowSOM result is updated to use this number of meta-clusters.
#' @param update_flowsom Logical. If `TRUE`, the `flowsom_result` will be updated in place with the new meta-clustering.
#' @param missing_seed Numeric.
#' If `flowsom_result` does not contain a "seed" element, this seed will be used for re-metaclustering.
#'
#' @return A named list with:
#' \describe{
#'   \item{`flowsom_newdata`}{The (optionally updated) `flowsom_result` with meta-clusters.}
#'   \item{`cells_clusters_from_train`}{Input `clustered_df` with meta-cluster annotations (if `clustered_df` was provided).}
#'   \item{`ncells_per_x`}{A summary `data.frame` (or list of `data.frame`s) with counts per meta-cluster (and cluster if applicable).}
#' }
#' @export
#' @keywords flowsom
#' ff_example <- example_processed()
#' fsom <- do_flowsom_TESTING(ff_example)
#' n_metaclusters <- 3
#' res_norecluster <- flowSOM_newMetacluster(
#'     flowsom_result = fsom,
#'     clustered_df = res[["ncells_per_x"]][["cluster"]],
#'     n_metacluster = n_metaclusters
#' )
flowSOM_newMetacluster <- function(flowsom_result,
                                   clustered_df = NULL,
                                   n_metacluster = NULL,
                                   update_flowsom = TRUE,
                                   missing_seed = 3711283) {
    # to avoid R CMD check note about undefined global variable
    cluster_numeric <- cluster <- value <- metaCluster <- ncells <- NULL
    if (!all(is.null(n_metacluster))) {
        # If the number of metaclusters is given, CHANGE the metaclusters to the new number
        if ("seed" %in% names(flowsom_result)) {
            seed <- flowsom_result$seed
        } else {
            seed <- missing_seed
        }

        if (n_metacluster <= 2) {
            stop("n_metacluster must be greater than 2, otherwise FlowSOM fails.")
        }

        if (!"ConsensusClusterPlus" %in% names(flowsom_result) ||
            length(flowsom_result[["ConsensusClusterPlus"]]) < n_metacluster) {
            # metaClustering_consensus() is a function from FlowSOM, only
            # calling ConsensusClusterPlus::ConsensusClusterPlus() and then returning the consensusClass
            # See the function here:
            # https://github.com/SofieVG/FlowSOM/blob/56892b583a0dbae5535665e45290ffce355bd503/R/4_metaClustering.R#L123
            # cl <- as.factor(FlowSOM::metaClustering_consensus(flowsom_result$map$codes,2 n_metacluster, seed = seed))
            cl_ccp <- suppressMessages(
                ConsensusClusterPlus::ConsensusClusterPlus(
                    t(flowsom_result$map$codes),
                    maxK = n_metacluster, reps = 100,
                    pItem = 0.9, pFeature = 1,
                    title = tempdir(), plot = "pdf", verbose = FALSE,
                    clusterAlg = "hc", distance = "euclidean", seed = seed
                )
            )

            # Store consensus class assignments in flowsom_result
            flowsom_result[["ConsensusClusterPlus"]] <- cl_ccp
            flowsom_result[["ConsensusClusterPlus_MAP"]] <- as.data.frame(sapply(cl_ccp[-1], function(x) {
                x$consensusClass
            }))
            colnames(flowsom_result[["ConsensusClusterPlus_MAP"]]) <- paste0("metaCluster_", 2:n_metacluster)
            flowsom_result[["ConsensusClusterPlus_MAP"]][, "cluster"] <- 1:nrow(flowsom_result[["ConsensusClusterPlus_MAP"]])
            flowsom_result[["ConsensusClusterPlus_MAP"]] <- flowsom_result[["ConsensusClusterPlus_MAP"]] |>
                dplyr::relocate(cluster)
        }

        # The as.factor comes from FlowSOM.
        # Assign meta-clusters (as factors) from stored mapping
        cl <- as.factor(flowsom_result[["ConsensusClusterPlus_MAP"]][[n_metacluster]])

        if (update_flowsom) {
            # Update FlowSOM object to reflect new meta-clustering
            flowsom_result$map$nMetaclusters <- length(levels(cl))
            flowsom_result$metaclustering <- cl
            # See what happens in the UpdateDerivedValues function:
            # https://github.com/SofieVG/FlowSOM/blob/56892b583a0dbae5535665e45290ffce355bd503/R/2_buildSOM.R#L116
            flowsom_result <- UpdateDerivedValues(flowsom_result)
            # Essentially, to assign metaClusters to (original) clusters, the following is done:
            #   flowsom_result$metaclustering[flowsom_result$map$mapping[, 1]]
            # So, flowsom_result$metaclustering is the metaCluster for each of the length(flowsom_result$metaclustering)
            # clusters.
        }
    }

    ncells_per_x <- NULL

    if (!is.null(clustered_df)) {
        # Either there must be a column called "cluster", then it is probably an already clustered (and previously exported)
        # data.frame of n cells and p channels, plus "cluster" and "metaCluster" columns.
        if ("cluster" %in% colnames(clustered_df)) {
            # Case 1: data.frame contains cluster assignments per cell
            clustered_df[["metaCluster"]] <- flowsom_result$metaclustering[clustered_df[["cluster"]]]

            ncells_per_x <- n_cells_per_x(
                dt = clustered_df,
                which = c("cluster", "metaCluster"),
                id_cols = c("sample"),
                which_levels = list(
                    "cluster" = as.character(1:flowsom_result$map$nNodes),
                    "metaCluster" = levels(flowsom_result$metaclustering)
                )
            )
        } else {
            # Case 2: clustered_df contains wide-format cluster columns (e.g., cluster_1, cluster_2, ...)
            clustered_df_long <- tidyr::pivot_longer(
                clustered_df,
                cols = tidyr::starts_with("cluster_"),
                names_to = "cluster",
                values_to = "value"
            ) |>
                dplyr::mutate(
                    cluster_numeric = as.numeric(stringr::str_remove(cluster, "cluster_"))
                )

            map_cluster_metacluster <- flowsom_result[["ConsensusClusterPlus_MAP"]][, c(1, n_metacluster)]
            colnames(map_cluster_metacluster) <- c("cluster", "metaCluster")

            clustered_df_long <- dplyr::left_join(
                clustered_df_long,
                map_cluster_metacluster,
                by = c("cluster_numeric" = "cluster")
            ) |>
                dplyr::select(-cluster_numeric) |>
                dplyr::summarise(
                    ncells = sum(value),
                    .by = c("sample", "metaCluster")
                ) |>
                dplyr::mutate(
                    metaCluster = paste0("metaCluster_", metaCluster)
                )

            metaclustered_df <- tidyr::pivot_wider(
                clustered_df_long,
                names_from = metaCluster,
                values_from = ncells,
                values_fill = 0
            )

            ncells_per_x <- list(
                "cluster" = clustered_df,
                "metaCluster" = clustered_df |>
                    dplyr::select(-tidyr::starts_with("cluster_")) |>
                    dplyr::left_join(metaclustered_df, by = "sample")
            )
        }
    }

    return(
        list(
            "flowsom_newdata" = flowsom_result, # Potentially updated to have the specified number of metaClusters
            "cells_clusters_from_train" = clustered_df,
            "ncells_per_x" = ncells_per_x
        )
    )
}
