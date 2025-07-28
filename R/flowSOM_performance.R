#' @title Evaluate Silhouette Performance from FlowSOM Clusterings
#'
#' @description
#' Computes silhouette scores to assess clustering performance using
#' FlowSOM results. Accepts both raw cluster assignments and the
#' result object from `flowSOM_optimal()`.
#'
#' @param dt_clustered data.table
#'   Either the clustered data table or a `flowSOM_optimal` result.
#' @param flowsom_result list
#'   The result from `flowSOM_optimal()`. Ignored if `dt_clustered`
#'   is already of class `flowSOM_optimal`.
#' @param n_metaclusters numeric (optional)
#'   Specifies which metacluster solution(s) to evaluate. Defaults
#'   to all available solutions in `flowsom_result`.
#' @param include_clusters logical
#'   Whether to include original cluster labels in the silhouette evaluation.
#' @param ncells integer
#'   Number of cells to subsample for silhouette score calculation.
#'   Defaults to 1000.
#' @param seed integer
#'   Random seed for reproducible subsampling.
#' @param relevant_cols character vector
#'   Column names used for distance calculation. If NULL, inferred from
#'   `flowsom_result`.
#'
#' @return data.table
#'   A data.table with silhouette scores for each clustering solution,
#'   including subsampled row indices.
flowSOM_performance <- function(
    dt_clustered,
    flowsom_result,
    n_metaclusters = NULL,
    include_clusters = TRUE,
    ncells = 1e3,
    seed = 42,
    relevant_cols = NULL) {
    # Support direct call: flowSOM_performance(flowsom_result)
    if ("flowSOM_optimal" %in% class(dt_clustered)) {
        flowsom_result <- dt_clustered
        dt_clustered <- flowsom_result[["cells_clusters_from_train"]]
    }
    dt_clustered <- data.table::data.table(dt_clustered)

    # Determine relevant columns for distance calculation
    if (is.null(relevant_cols)) {
        if (missing(flowsom_result)) {
            # Heuristically exclude known metadata columns
            exclude_cols <- c("sample", "is_som_outlier", "cluster", "metaCluster", "Time")
            relevant_cols <- setdiff(colnames(dt_clustered), exclude_cols)
            warning("No relevant_cols provided. Using all columns except known metadata columns.")
        } else {
            relevant_cols <- flowsom_result$fs_res_train$map$colsUsed
        }
    }

    # Determine metacluster mappings to evaluate
    n_metaclusters <- tryCatch(
        seq_len(flowsom_result$fs_res_train$map$nMetaclusters),
        error = function(e) NULL
    )

    # Prepare cluster-to-metacluster mapping
    if (!is.null(n_metaclusters)) {
        clustering_map <- flowsom_result$fs_res_train$ConsensusClusterPlus_MAP[, n_metaclusters, drop = FALSE]
    } else {
        # Fall back to unique cluster–metaCluster pairs
        clustering_map <- unique(dt_clustered[, .(cluster, metaCluster)])[, metaCluster := as.numeric(as.character(metaCluster))]
    }

    # Optional subsampling of cells for efficiency
    dt_clustered[, rownum := seq_len(.N)]
    if (!is.null(ncells)) {
        dt_clustered[, rownum := seq_len(.N)]
        set.seed(seed)
        dt_clustered <- dt_clustered[sample(.N, min(ncells, .N))]
    }

    # Select appropriate mapping(s) to evaluate
    test_clustermap <- if (include_clusters) clustering_map else clustering_map[, -1, with = FALSE]

    # Compute distance matrix from selected features
    distmat <- distances::distances(
        dt_clustered[, ..relevant_cols]
    )

    # Apply silhouette calculation for each clustering solution
    counter <- 0
    score_sil <- apply(test_clustermap, 2, function(clusters) {
        counter <<- counter + 1
        # Remap cluster to metacluster using current mapping
        current_map <- data.table(cluster = clustering_map$cluster, metaCluster = clusters)

        # Join metaCluster labels to dt_clustered
        meta_labels <- current_map[dt_clustered[, .(cluster)], on = "cluster"][["metaCluster"]]
        # Compute silhouette scores
        sil <- cluster::silhouette(meta_labels, dist = distmat) |>
            data.table::as.data.table()
        colnames(sil) <- c("metacluster", "neighbor", "sil_width")
        # Return silhouette results along with subsampled row indices
        cbind(rownum = dt_clustered[["rownum"]], sil)
    }) |> data.table::rbindlist(idcol = "clustering")

    return(score_sil)
}
