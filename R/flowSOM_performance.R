#' @title Get performances from flowSOM clusterings
#' @param flowsom_result
#'  Result of flowSOM_optimal()
flowSOM_performance <- function(
    dt_clustered,
    flowsom_result,
    n_metaclusters = NULL,
    include_clusters = TRUE,
    ncells = 1e3,
    seed = 42,
    relevant_cols = NULL) {
    if ("flowSOM_optimal" %in% class(dt_clustered)) {
        # then the dt_clustered is the result of flowSOM_optimal()
        # and we should regard dt_clustered as flowsom_result.
        # This is such that flowSOM_performance(dt_clustered) but also
        # flowSOM_performance(flowsom_result) works.
        flowsom_result <- dt_clustered
        dt_clustered <- flowsom_result[["cells_clusters_from_train"]]
    }

    if (is.null(relevant_cols)) {
        if (missing(flowsom_result)) {
            relevant_cols <- colnames(dt_clustered)[!colnames(dt_clustered) %in% c("sample", "is_som_outlier", "cluster", "metaCluster", "Time")]
            warning("No relevant_cols provided, using all columns except 'sample', 'is_som_outlier', 'cluster', 'metaCluster', and 'Time'.")
        } else {
            relevant_cols <- flowsom_result$fs_res_train$map$colsUsed
        }
    }
    # If the following doesnt work, then it is still NULL (and handled right after.)
    n_metaclusters <- tryCatch(
        seq_len(flowsom_result$fs_res_train$map$nMetaclusters),
        error = function(e) {
            NULL
        }
    )

    clustering_map <- NULL
    if (!all(is.null(n_metaclusters))) {
        clustering_map <- flowsom_result$fs_res_train$ConsensusClusterPlus_MAP[, n_metaclusters, drop = FALSE]
    } else {
        # clustering map should only use the existing metaCluster column
        clustering_map <- dt_clustered[, c("cluster", "metaCluster"), with = FALSE] |> unique()
    }

    if (!is.null(ncells)) {
        dt_clustered[, rownum := seq_len(.N)]
        set.seed(seed)
        dt_clustered <- dt_clustered[sample(.N, min(ncells, .N)), ]
    }


    if (!include_clusters) {
        test_clustermap <- clustering_map[, -1]
    } else {
        test_clustermap <- clustering_map
    }
    distmat <- distances::distances(
        dt_clustered[, relevant_cols, with = FALSE]
    )
    score_sil <- apply(
        test_clustermap, 2, function(clusters) {
            current_map <- data.table::data.table("cluster" = clustering_map[["cluster"]], "metaCluster" = clusters)

            sil <- cluster::silhouette(
                dt_clustered[, "cluster", with = FALSE][current_map, on = "cluster"][["metaCluster"]] |>
                    as.character() |> # in case it's a factor, first needs to be converted to character
                    as.numeric(),
                dist = distmat
            ) |>
                data.table::as.data.table()
            # If rownum does not exist, then the "new" rownum column is just not present thanks to NULL
            return(
                cbind("rownum" = dt_clustered[["rownum"]], sil)
            )
        }
    ) |> data.table::rbindlist(idcol = "clustering")

    return(score_sil)
}
