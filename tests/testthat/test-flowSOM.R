test_that("FlowSOM outliers", {
    ff_example <- example_processed(ncells = Inf)
    fsom <- do_flowsom_TESTING(ff_example)
    outliers <- flowSOM_is.outlier(fsom)
    testthat::expect_equal(sum(outliers), 406)
    testthat::expect_length(outliers, flowCore::nrow(ff_example))
})


test_that("FlowSOM outliers 2", {
    ff_example <- example_processed(ncells = Inf)
    fsom <- do_flowsom_TESTING(ff_example)
    res <- flowSOM_predict(flowsom_result = fsom, flowset = flowCore::flowSet(ff_example))

    outliers <- flowSOM_is.outlier(fsom)
    testthat::expect_equal(outliers, res$cells_clusters_from_train[["is_som_outlier"]])
    testthat::expect_equal(sum(outliers), 406)
    testthat::expect_equal(sum(res$cells_clusters_from_train[["is_som_outlier"]]), 406)
})

test_that("FlowSOM new Metaclusters", {
    ff_example <- example_processed()
    fsom <- do_flowsom_TESTING(ff_example)
    testthat::expect_true(fsom$map$nMetaclusters == 10)

    # for (n_metaclusters in 3:44) {
    for (n_metaclusters in c(3, 8, 12)) {
        res <- flowSOM_predict(
            flowsom_result = fsom,
            flowset = flowCore::flowSet(ff_example),
            n_metacluster = n_metaclusters
        )
        testthat::expect_true(res$flowsom_newdata$map$nMetaclusters == n_metaclusters)
        testthat::expect_true(ncol(res[["ncells_per_x"]][["metaCluster"]]) == n_metaclusters + 1)
    }
    testthat::expect_error(
        flowSOM_predict(
            flowsom_result = fsom,
            flowset = flowCore::flowSet(ff_example),
            n_metacluster = 2
        )
    )
    testthat::expect_error(
        flowSOM_predict(
            flowsom_result = fsom,
            flowset = flowCore::flowSet(ff_example),
            n_metacluster = 1
        )
    )
    testthat::expect_error(
        flowSOM_predict(
            flowsom_result = fsom,
            flowset = flowCore::flowSet(ff_example),
            n_metacluster = -1
        )
    )
    testthat::expect_error(
        # maximum seems to be 45 in this case from flowSOM.
        flowSOM_predict(
            flowsom_result = fsom,
            flowset = flowCore::flowSet(ff_example),
            n_metacluster = 45
        )
    )
})
# test wrapper_count_model
test_that("FlowSOM wrapper_count_model", {
    ff_example <- example_processed()
    fsom <- do_flowsom_TESTING(ff_example)
    res <- flowSOM_predict(
        flowsom_result = fsom,
        flowset = flowCore::flowSet(ff_example)
    )
    fake_clusterings <- lapply(res[["ncells_per_x"]], function(x) {
        tmp <- lapply(1:25, function(y) {
            x[1, -1] <- as.list(sample(1000, ncol(x) - 1))
            return(x)
        }) |> do.call(what = rbind)
        set.seed(40)
        tmp[["tvt"]] <- sample(
            c("train", "validation", "test", "prospective"),
            nrow(tmp),
            replace = TRUE
        )
        tmp[["sample"]] <- sample(c("A", "B"), nrow(tmp), replace = TRUE)
        return(tmp)
    })
    allowed_clusterings <- c("cluster", "metaCluster")
    w <- testthat::capture_warnings(
        finalmodels_predictions <- wrapper_count_models(
            df_list = fake_clusterings[allowed_clusterings],
            tvt_col = "tvt",
            # outdir = file.path(outdir_base, basename(metaclustering_dir), training_device),
            outdir = local_tempdir_time(),
            dvs_potential = c("sample"),
            dvs_multiclass = c(),
            ivs_regex = "[cC]luster",
            hparam_n_evaluations = 3,
            seed = 1372873,
            learners_classification = list(
                mlr3::lrn(
                    "classif.ranger",
                    predict_type = "prob", predict_sets = c("train", "test"),
                    max.depth = paradox::to_tune(2, 20), # minimum and maximum depth
                    num.trees = paradox::to_tune(c(500, 1000, 1500, 2000)),
                    importance = "impurity"
                )
            ),
            dv_class_positive = c("sample" = "B"),
            measures = mlr3::msr("classif.logloss")
        )
    )
    testthat::expect_match(w[1], "already exists", all = TRUE)
    # # I think the following is not on me, therefore I will not fix it
    # testthat::expect_match(w[-1], "weights' is deprecated.*Use.*weights_learner", all = TRUE)
})

test_that("FlowSOM new Metaclusters, with/out new clustering", {
    ff_example <- example_processed()
    fsom <- do_flowsom_TESTING(ff_example)
    testthat::expect_true(fsom$map$nMetaclusters == 10)

    for (n_metaclusters in c(3, 8, 12)) {
        res <- flowSOM_predict(
            flowsom_result = fsom,
            flowset = flowCore::flowSet(ff_example),
            n_metacluster = n_metaclusters
        )

        res_norecluster <- flowSOM_newMetacluster(
            flowsom_result = fsom,
            clustered_df = res[["ncells_per_x"]][["cluster"]],
            n_metacluster = n_metaclusters
        )
        testthat::expect_true(res$flowsom_newdata$map$nMetaclusters == n_metaclusters)
        testthat::expect_true(ncol(res[["ncells_per_x"]][["metaCluster"]]) == n_metaclusters + 1)
        testthat::expect_equal(res["ncells_per_x"], res_norecluster["ncells_per_x"])
        testthat::expect_false(identical(res_norecluster[["cells_clusters_from_train"]], res[["cells_clusters_from_train"]]))


        res_norecluster_cellwise <- flowSOM_newMetacluster(
            flowsom_result = fsom,
            clustered_df = res[["cells_clusters_from_train"]],
            n_metacluster = n_metaclusters
        )
        testthat::expect_equal(res["ncells_per_x"], res_norecluster_cellwise["ncells_per_x"])
        testthat::expect_true(identical(res_norecluster_cellwise[["cells_clusters_from_train"]], res[["cells_clusters_from_train"]]))
        for (relevant_value in c(
            "pattern", "compensate", "spillover", "transform", "toTransform",
            "transformFunction", "transformList", "scale", "metaData", "map",
            "MST", "metaclustering", "outliers"
        )) {
            testthat::expect_identical(res_norecluster_cellwise[[relevant_value]], res[[relevant_value]])
            testthat::expect_identical(res_norecluster[[relevant_value]], res[[relevant_value]])
        }
    }
})

test_that("FlowSOM optimal", {
    ff_example <- example_processed()
    testthat::expect_warning(
        fsom <- flowSOM_optimal(
            flowCore::flowSet(ff_example),
            # Input options:
            compensate = FALSE,
            transform = FALSE,
            scale = FALSE,
            # SOM options:
            colsToUse = c(9, 12, 14:18), xdim = 7, ydim = 7,
            # Metaclustering options:
            nClus = 10
        ),
        "By using nClus, I am ignoring the parameter maxMeta"
    )

    testthat::expect_true(fsom$fs_res_train$map$nMetaclusters == 10)
    testthat::expect_equal(dim(fsom$fs_res_train$ConsensusClusterPlus_MAP), c(49, 10))
})

test_that("FlowSOM optimal, defaults", {
    ff_example <- example_processed()
    testthat::expect_warning(
        suppressMessages(
            fsom <- flowSOM_optimal(
                flowCore::flowSet(ff_example),
                # Metaclustering options:
                nClus = 10
            )
        ),
        "By using nClus, I am ignoring the parameter maxMeta"
    )

    testthat::expect_true(fsom$fs_res_train$map$nMetaclusters == 10)
    testthat::expect_equal(dim(fsom$fs_res_train$ConsensusClusterPlus_MAP), c(100, 10))
})

test_that("FlowSOM optimal, outdir", {
    ff_example <- example_processed()
    tmpdir <- local_tempdir_time()
    testthat::expect_warning(
        suppressMessages(
            fsom <- flowSOM_optimal(
                flowCore::flowSet(ff_example),
                # Metaclustering options:
                nClus = 10,
                outdir = tmpdir
            )
        ),
        "By using nClus, I am ignoring the parameter maxMeta"
    )
    loaded_fsomresult2 <- qs2::qs_read(file.path(tmpdir, "r1-FlowSOM_result_train.qs2"))

    testthat::expect_equivalent(loaded_fsomresult2, fsom$fs_res_train)
    testthat::expect_true(fsom$fs_res_train$map$nMetaclusters == 10)
    testthat::expect_equal(dim(fsom$fs_res_train$ConsensusClusterPlus_MAP), c(100, 10))
})
