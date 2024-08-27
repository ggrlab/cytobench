#' @title Predict new cells using a flowsom result
#' @description
#' Predict the (meta-) clusters for all cells using a flowsom result.
#' @param flowsom_result Result of flowSOM_optimal() or flowSOM()
#' @param flowset The flowset whose cells should be assigned to the clusters from flowsom_result
#' @param madAllowed
#' See FlowSOM::TestOutliers() or flowSOM_is.outlier()
#' @export
flowSOM_predict <- function(flowsom_result, flowset, madAllowed = 4) {
    if ("fs_res_train" %in% names(flowsom_result)) {
        flowsom_result <- flowsom_result[["fs_res_train"]]
    }
    colsToUse <- flowsom_result$map$colsUsed
    ### 3. Predict the (meta-) clusters for all cells (training AND test set)
    # 3.1 Bind all cells from all samples
    fcs_first_dt_long <- flowCore::fsApply(flowset, function(y) {
        data.table::as.data.table(flowCore::exprs(y))
    }, simplify = FALSE) |>
        data.table::rbindlist(idcol = "sample", fill = FALSE)

    # 3.2 Predict the clusters for all cells
    predicted_fs_train_allcells <- FlowSOM::NewData(
        flowsom_result,
        input = as.matrix(fcs_first_dt_long[, colsToUse, with = FALSE]),
        # If transform and scale are NULL, the same as in the flowSOM function will be used
    )

    cells_clusters_from_train <- data.table::data.table(fcs_first_dt_long)
    cells_clusters_from_train[, metaCluster := FlowSOM::GetMetaclusters(predicted_fs_train_allcells)]
    cells_clusters_from_train[, cluster := FlowSOM::GetClusters(predicted_fs_train_allcells)]

    all_cluster_x_ids <- list(
        "cluster" = as.character(1:predicted_fs_train_allcells$map$nNodes),
        "metaCluster" = levels(predicted_fs_train_allcells$metaclustering)
    )
    # # 3.3 Save the predicted clusters
    # if (!is.null(outdir)) {
    #     qs::qsave(cells_clusters_from_train, file.path(outdir, "r2-FlowSOM_predicted_withTrain.qs"))
    # }
    ### 4. Extract the number of cells per cluster
    # id_cols <- c("tvt", "sample")
    id_cols <- c("sample")
    ncells_per_x <- sapply(c("cluster", "metaCluster"), function(x) {
        grouping_columns <- c(id_cols, x)
        tmp <- cells_clusters_from_train[, .N, by = grouping_columns]
        tmp[[x]] <- factor(tmp[[x]], levels = sort(as.numeric(all_cluster_x_ids[[x]])))
        levels(tmp[[x]]) <- paste0(x, "_", levels(tmp[[x]]))

        removin_sample <- data.frame(
            "_____REMOVEME_____", levels(tmp[[x]]), 0
        )
        # ensure that the column names are the same
        colnames(removin_sample) <- colnames(tmp)
        tmp_wide <- dplyr::add_row(
            tmp,
            # Add a REMOVEME sample just to make sure that all clusters are present in the output
            removin_sample
        ) |>
            tidyr::pivot_wider(
                names_from = tidyr::all_of(x),
                values_from = "N",
                values_fill = 0
            ) |>
            # After the pivot_wider, the REMOVEME sample is not needed anymore
            dplyr::filter(sample != "_____REMOVEME_____")
        # The following mainly resorts the cluster_I columns
        tmp_wide[, c(id_cols, levels(tmp[[x]]))]
    }, USE.NAMES = TRUE, simplify = FALSE)

    ### 5. Get the outliers
    is_cell_outlier <- flowSOM_is.outlier(
        fsom = predicted_fs_train_allcells,
        fsomReference = flowsom_result,
        madAllowed = madAllowed
    )
    cells_clusters_from_train[["is_som_outlier"]] <- is_cell_outlier
    res <- list(
        cells_clusters_from_train = cells_clusters_from_train,
        ncells_per_x = ncells_per_x,
        flowsom_newdata = predicted_fs_train_allcells
    )
    return(res)
}
