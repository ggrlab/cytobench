#' @title FlowSOM Clustering with Optimal Number of MetaClusters
#'
#' @description
#' Runs the FlowSOM algorithm on a cytometry dataset and determines an optimal number
#' of meta-clusters using `ConsensusClusterPlus`. Automatically generates and optionally
#' saves plots and the clustering result. Can also return per-cell cluster predictions.
#'
#' @param fs_train A `flowSet` or `cytoset` representing the cytometry training data.
#' @param colsToUse Character vector. Column names (channels/markers) to include in clustering.
#' @param transformList A `transformList` (typically from `flowCore::transformList`) defining the transformations to apply.
#'        If `transform = TRUE`, this will be applied during clustering.E.g.
#'
#' \code{
#' fc_transforms <- sapply(colsToUse, function(cn_x) {
#'     flowCore::arcsinhTransform(
#'         a = 0,
#'         b = 1 / pre_rescale_asinh[["navios"]][[cn_x]],
#'         c = 0
#'     )
#' }, simplify = FALSE)
#' fc_transformlist_asinh <- flowCore::transformList(names(fc_transforms), fc_transforms)
#' }
#' @param outdir Optional. Directory path to save the FlowSOM result object and associated PDF plots.
#' @param seed Integer. Random seed for reproducibility. Default is `3711283`.
#' @param scale Logical. Whether to scale data (default `TRUE`). Scaling is applied prior to SOM training.
#' @param transform Logical. Whether to apply `transformList` before clustering (default `TRUE`).
#' @param maxMeta Integer. Maximum number of meta-clusters to consider. If `NULL`, defaults to a FIXED number via `nClus` argument.
#' @param ... Additional arguments passed to `FlowSOM::FlowSOM()`.
#'
#' @return An object of class `"flowSOM_optimal"`, which includes:
#' \describe{
#'   \item{`fs_res_train`}{The FlowSOM result object with meta-clustering annotations.}
#'   \item{...}{Additional elements returned by `flowSOM_predict()`, including predicted cluster and meta-cluster labels.}
#' }
#'
#' @export
#' @keywords flowsom
#' @examples
#' ff_example <- example_processed()
#' fsom <- flowSOM_optimal(
#'     flowCore::flowSet(ff_example),
#'     # Input options:
#'     compensate = FALSE,
#'     transform = FALSE,
#'     scale = FALSE,
#'     # SOM options:
#'     colsToUse = c(9, 12, 14:18), xdim = 3, ydim = 3,
#'     # Metaclustering options:
#'     nClus = 5
#' )
flowSOM_optimal <- function(fs_train,
                            colsToUse = NULL,
                            transform = FALSE,
                            transformList = NULL,
                            outdir = NULL,
                            seed = 3711283,
                            scale = TRUE,
                            maxMeta = 50, # Max number of metaClusters
                            ...) {
    # Warn if very few meta-clusters are requested
    if (is.null(maxMeta) || maxMeta < 10) {
        warning("Are you sure you want to set maxMeta to ", maxMeta, "? This is the maximal number of possible metaclusters.")
    }

    # If user provides nClus directly, override maxMeta and show warning
    if ("nClus" %in% names(list(...))) {
        warning("By using nClus, I am ignoring the parameter maxMeta!")
        maxMeta <- NULL
    }

    # 1. Run the FlowSOM clustering algorithm
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

    # 2. Precompute all meta-clusterings for ConsensusClusterPlus
    # With this, we can later very quickly get all meta-clusterings up until n_meta using
    # flowSOM_newMetacluster(fs_res_train_allmeta, n_metacluster = XXXX)
    n_meta <- max(maxMeta, list(...)[["nClus"]])
    fs_res_train_allmeta <- flowSOM_newMetacluster(
        fs_res_train,
        n_metacluster = n_meta,
        update_flowsom = FALSE
    )
    fs_res_train <- fs_res_train_allmeta[["flowsom_newdata"]]

    # 3. Save result and FlowSOM plots if requested
    if (!is.null(outdir)) {
        dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
        qs2::qs_save(fs_res_train, file.path(outdir, "r1-FlowSOM_result_train.qs2"))

        grDevices::pdf(file.path(outdir, "p1-FlowSOM.pdf"), width = 10, height = 10)
        FlowSOM::PlotStars(
            fs_res_train,
            backgroundValues = fs_res_train$metaclustering
        )
        FlowSOM::PlotNumbers(fs_res_train)
        FlowSOM::PlotNumbers(fs_res_train, level = "metaclusters")

        for (marker_x in colsToUse) {
            print(
                FlowSOM::PlotMarker(
                    fs_res_train,
                    marker = marker_x,
                    backgroundValues = fs_res_train$metaclustering
                )
            )
        }
        grDevices::dev.off()
    }

    # 4. Predict the per-cell cluster and metacluster assignment
    predicted_fs <- flowSOM_predict(
        flowsom_result = fs_res_train,
        flowset = fs_train
    )

    # 5. Combine and return results
    res <- c(
        list(fs_res_train = fs_res_train),
        predicted_fs
    )
    class(res) <- c("flowSOM_optimal")
    return(res)
}
