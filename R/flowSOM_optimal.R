#' @title flowSOM clustering with optimal number of clusters
#' @description
#' Clusters the cells using the FlowSOM algorithm with the optimal number of clusters (determined by ConsensusClusterPlus), by setting maxMeta.
#' @param fs_train
#' The training set of the cytometry data.
#' @param colsToUse
#' The columns to use for clustering.
#' @param transformList
#' The list of transformations to apply to the data. E.g.
#'
#' fc_transforms <- sapply(colsToUse, function(cn_x) {
#'     flowCore::arcsinhTransform(
#'         a = 0,
#'         b = 1 / pre_rescale_asinh[["navios"]][[cn_x]],
#'         c = 0
#'     )
#' }, simplify = FALSE)
#' fc_transformlist_asinh <- flowCore::transformList(names(fc_transforms), fc_transforms)
#' @param outdir
#' The directory to save the FlowSOM results.
#' @param seed
#' The seed for the random number generator.
#' @param scale
#' Whether to scale the data.
#' @param transform
#' Whether to transform the data using transformList
#' @param maxMeta
#' The maximal number of metaclusters. If NULL, the number of clusters is FIXED by nClus.
#' @param ... Further parameters to FlowSOM
#' @export
flowSOM_optimal <- function(fs_train,
                            colsToUse,
                            transform = TRUE,
                            transformList = NULL,
                            outdir = NULL,
                            seed = 3711283,
                            scale = TRUE,
                            maxMeta = 50, # Number of MetaClusters, automatically determined by ConsensusClusterPlus
                            ...) {
    if (is.null(maxMeta) || maxMeta < 10) {
        warning("Are you sure you want to set maxMeta to ", maxMeta, "? This is the maximal number of possible metaclusters.")
    }
    if ("nClus" %in% names(list(...))) {
        warning("By using nClus, I am ignoring the parameter maxMeta!")
        maxMeta <- NULL
    }
    ### 1. Run FlowSOM
    fs_res_train <- FlowSOM::FlowSOM(
        input = fs_train,
        transform = transform,
        transformList = transformList,
        colsToUse = colsToUse,
        maxMeta = maxMeta,
        scale = scale,
        seed = seed,
        ...
    )
    fs_res_train[["seed"]] <- seed
    # The following creates the metacluster dataframe for ALL metaclusterings
    #  fs_res_train$ConsensusClusterPlus_MAP
    fs_res_train_allmeta <- flowSOM_newMetacluster(fs_res_train, n_metacluster = max(maxMeta, list(...)[["nClus"]]), update_flowsom = FALSE)
    fs_res_train <- fs_res_train_allmeta[["flowsom_newdata"]]

    ### 2. Plot and save FlowSOM results into outdir
    if (!is.null(outdir)) {
        dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
        qs::qsave(fs_res_train, file.path(outdir, "r1-FlowSOM_result_train.qs"))

        pdf(file.path(outdir, "p1-FlowSOM.pdf"), width = 10, height = 10)
        FlowSOM::PlotStars(
            fs_res_train,
            backgroundValues = fs_res_train$metaclustering
        )
        FlowSOM::PlotNumbers(
            fs_res_train
        )
        FlowSOM::PlotNumbers(
            fs_res_train,
            level = "metaclusters"
        )
        for (marker_x in colsToUse) {
            print(
                FlowSOM::PlotMarker(
                    fs_res_train,
                    marker = marker_x,
                    backgroundValues = fs_res_train$metaclustering
                )
            )
        }
        dev.off()
    }

    predicted_fs <- flowSOM_predict(
        flowsom_result = fs_res_train,
        flowset = fs_train
    )

    return(c(
        list(
            fs_res_train = fs_res_train
        ), predicted_fs
    ))
}
