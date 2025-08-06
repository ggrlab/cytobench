test_that("FlowSOM optimal", {
    ff_example <- example_processed()
    fsom <- do_flowsom_TESTING(ff_example)
    w <- testthat::capture_warnings({
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
    })
    testthat::expect_match(w, "By using nClus, I am ignoring the parameter maxMeta", all = FALSE)
    testthat::expect_warning(
        scores1 <- flowSOM_performance(fsom_opt$cells_clusters_from_train), 
        "No relevant_cols provided. Using all columns except known metadata columns."
    )
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
