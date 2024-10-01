#' @title Predict new cells using a flowsom result
#' @description
#' Predict the (meta-) clusters for all cells using a flowsom result.
#' @param flowsom_result Result of flowSOM_optimal() or flowSOM()
#' @param clustered_df
#'  1) The data.frame including the clustered cells. The column "cluster" must be present and
#'  should be the cluster assignment from the FlowSOM ("flowsom_result").
#'
#' 2) Alternatively: The data.frame including the number of clustered cells. The columns "cluster_1", "cluster_2", ..., "cluster_n".
#' The function will then calculate the metaClusters for each of the clusters and return the number of cells per metaCluster.
#' @param n_metacluster
#' If not NULL, the number of metaclusters will be changed to this number in flowsom_result.
#' @return
#' A list with two elements:
#' 1) "cells_clusters_from_train": The data.frame with the new metaClusters. If the input was a data.frame with the number of cells per cluster,
#' the output will be a data.frame with the number of cells per metaCluster. (Identical to 2))
#' 2) "ncells_per_x": A data.frame in wide format with the counts of cells per specified grouping.
#' @export
flowSOM_newMetacluster <- function(flowsom_result,
                                   clustered_df,
                                   n_metacluster = NULL,
                                   update_flowsom = TRUE) {
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
        if (!"ConsensusClusterPlus" %in% names(flowsom_result) || length(flowsom_result[["ConsensusClusterPlus"]]) < n_metacluster) {

            # metaClustering_consensus() is a function from FlowSOM, only
            # calling ConsensusClusterPlus::ConsensusClusterPlus() and then returning the consensusClass
            # See the function here:
            # https://github.com/SofieVG/FlowSOM/blob/56892b583a0dbae5535665e45290ffce355bd503/R/4_metaClustering.R#L123
            # cl <- as.factor(FlowSOM::metaClustering_consensus(flowsom_result$map$codes,2 n_metacluster, seed = seed))
            cl_ccp <- suppressMessages(
                ConsensusClusterPlus::ConsensusClusterPlus(t(flowsom_result$map$codes),
                    maxK = n_metacluster, reps = 100, pItem = 0.9, pFeature = 1, title = tempdir(),
                    plot = "pdf", verbose = FALSE, clusterAlg = "hc", distance = "euclidean",
                    seed = seed
                )
            )
            flowsom_result[["ConsensusClusterPlus"]] <- cl_ccp
        }
        # The as.factor comes from FlowSOM.
        cl <- as.factor(flowsom_result[["ConsensusClusterPlus"]][[n_metacluster]]$consensusClass)

        if (update_flowsom) {
            # Update the flowSOM object
            flowsom_result$map$nMetaclusters <- length(levels(cl))
            flowsom_result$metaclustering <- cl
            # See what happens in the FlowSOM:::UpdateDerivedValues function:
            # https://github.com/SofieVG/FlowSOM/blob/56892b583a0dbae5535665e45290ffce355bd503/R/2_buildSOM.R#L116
            flowsom_result <- FlowSOM:::UpdateDerivedValues(flowsom_result)
            # Essentially, to assign metaClusters to (original) clusters, the following is done:
            #   flowsom_result$metaclustering[flowsom_result$map$mapping[, 1]]
            # So, flowsom_result$metaclustering is the metaCluster for each of the length(flowsom_result$metaclustering)
            # clusters.
        }
    }

    # Either there must be a column called "cluster", then it is probably an already clustered (and previously exported)
    # data.frame of n cells and p channels, plus "cluster" and "metaCluster" columns.
    if ("cluster" %in% colnames(clustered_df)) {
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
        # Then the clusters must be in the columns of the data.frame like "cluster_1", "cluster_2", ..., "cluster_n"
        clustered_df_long <- tidyr::pivot_longer(
            clustered_df,
            cols = tidyr::starts_with("cluster_"),
            names_to = "cluster",
            values_to = "value"
        ) |>
            dplyr::mutate(
                cluster_numeric = as.numeric(stringr::str_remove(cluster, "cluster_"))
            )

        map_cluster_metacluster <- data.frame(
            cluster = 1:flowsom_result$map$nNodes,
            metaCluster = flowsom_result$metaclustering[1:flowsom_result$map$nNodes]
        )

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
                dplyr::left_join(
                    metaclustered_df,
                    by = "sample"
                )
        )
    }

    return(
        list(
            "flowsom_newdata" = flowsom_result, # Potentially updated to have the specified number of metaClusters
            "cells_clusters_from_train" = clustered_df,
            "ncells_per_x" = ncells_per_x
        )
    )
}
