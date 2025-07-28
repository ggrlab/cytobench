example_processed <- function() {
    fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
    ff <- flowCore::read.FCS(fileName)
    # Compensation
    comp <- flowCore::keyword(ff)[["SPILL"]]
    ff <- flowWorkspace::compensate(ff, comp)
    # Transformation
    transformList <- flowCore::estimateLogicle(ff, channels = colnames(comp))
    ff <- flowWorkspace::transform(ff, transformList)
    return(ff)
}
do_flowsom <- function(ff) {
    set.seed(237123)
    fSOM <- FlowSOM::FlowSOM(ff,
        # Input options:
        compensate = FALSE,
        transform = FALSE,
        scale = FALSE,
        # SOM options:
        colsToUse = c(9, 12, 14:18), xdim = 7, ydim = 7,
        # Metaclustering options:
        nClus = 10
    )
    fSOM[["seed"]] <- 237123
    return(fSOM)
}

test_that("FlowSOM optimal", {
    ff_example <- example_processed()[1:1000, ]
    fsom <- do_flowsom(ff_example)
    fsom_opt <- flowSOM_optimal(
        flowCore::flowSet(ff_example),
        # Input options:
        compensate = FALSE,
        transform = FALSE,
        scale = FALSE,
        # SOM options:
        colsToUse = c(9, 12, 14:18), xdim = 3, ydim = 3,
        # Metaclustering options:
        nClus = 3
    )
    scores1 <- flowSOM_performance(fsom_opt$cells_clusters_from_train)
    scores1_v2 <- flowSOM_performance(
        fsom_opt$cells_clusters_from_train,
        relevant_cols = fsom_opt$fs_res_train$map$colsUsed,
    )
    scores2 <- flowSOM_performance(fsom_opt)
    scores3 <- flowSOM_performance(
        dt_clustered = fsom_opt$cells_clusters_from_train,
        flowsom_result = fsom_opt,
    )
    testthat::expect_true(identical(scores3, scores2))
    # Scores 1 must be different, because it uses the default relevant_cols (being all columns except known metadata columns)
    cols_to_be_equal <- c("rownum", "metacluster", "neighbor", "sil_width")
    testthat::expect_false(identical(
        scores1[, c(cols_to_be_equal), with = FALSE],
        scores2 |> dplyr::filter(clustering %in% c("cluster", "metaCluster_3")) |> dplyr::select(tidyr::all_of(cols_to_be_equal))
    ))
    testthat::expect_equal(
        scores1_v2[, c(cols_to_be_equal), with = FALSE],
        scores2 |> dplyr::filter(clustering %in% c("cluster", "metaCluster_3")) |> dplyr::select(tidyr::all_of(cols_to_be_equal))
    )
})
