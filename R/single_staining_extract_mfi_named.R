#' Extract Median Fluorescence Intensities from single-stained samples
#'
#' Extract 
#'
#' @param single_stainings DESCRIPTION.
#' @param sample_unstained DESCRIPTION.
#' @param relevant_difference DESCRIPTION.
#' @param minimum_fraction_per_cluster DESCRIPTION.
#'
#' @return RETURN_DESCRIPTION
#' @examples
#' # ADD_EXAMPLES_HERE
#' 
single_staining_extract_mfi_named <- function(
    single_stainings,
    sample_unstained = NULL,
    ignore_features = c("FS INT", "FS TOF", "SS INT", "TIME"),
    minimum_diff = 5e3) {
    # single_stainings:
    #   Named list of matrices, each including one single-stained value
    #   The names of the list should be the names of the singly-stained features

    ss_clusterings <- lapply(names(single_stainings), function(feature) {
        mat <- single_stainings[[feature]]
        if (!feature %in% colnames(mat)) {
            stop(paste0("Feature '", feature, "' not found in the matrix."))
        }

        ss_feature_fit <- kmeans(mat[[feature]], centers = 2)

        ## Testing various methods of clustering to see if that makes a difference
        # nolint start
        # pdf("removeme_2.pdf")
        # ggplot(data.frame(val = mat[[feature]], cluster = ss_feature_fit$cluster)) +
        #     geom_histogram(aes(x = val, fill = factor(cluster)), alpha=.6, bins=50, position="identity")
        # dev.off()

        # browser()
        # ss_feature_fit <- fastcluster::hclust(dist(mat[[feature]]))
        # ss_feature_fit_cut <- cutree(ss_feature_fit, k=2)
        # pdf("removeme_2.pdf")
        # ggplot(data.frame(val = mat[[feature]], cluster = ss_feature_fit_cut)) +
        #     geom_histogram(aes(x = val, fill = factor(cluster)), alpha=.6, bins=50, position="identity")
        # dev.off()
        # pacman::p_load(dbscan)
        # # tmp <- dbscan::hdbscan(as.matrix(mat[[feature]]), minPts = 3000)
        # tmp <- dbscan::dbscan(as.matrix(mat[[feature]]), eps = .1)
        # pdf("removeme_2.pdf")
        # ggplot(data.frame(val = mat[[feature]], cluster = tmp$cluster)) +
        #     geom_histogram(aes(x = val, fill = factor(cluster)), alpha=.6, bins=50, position="identity") +
        #     theme(legend.position = "none")
        # dev.off()

        # pacman::p_load("mclust")
        # # Clustering
        # mod1 <- mclust::Mclust(mat[[feature]])
        # summary(mod1)
        # plot(mod1,  what = c("BIC", "classification"))
        # nolint end

        cluster_medians <- aggregate(
            mat,
            by = list(ss_feature_fit$cluster),
            FUN = median
        )
        if (abs(diff(cluster_medians[, feature])) < minimum_diff) {
            stop(
                "The difference between the two clusters according to the feature '",
                feature, "' is lower than ", minimum_diff,
                "(", abs(diff(cluster_medians[, feature])),
                "). Is this a single-stained sample in linear(per default settings) space?"
            )
        }
        # remove the introduced cluster column "Group.1"
        retmat <- as.matrix(sort(cluster_medians[, feature]))
        colnames(retmat) <- feature
        return(retmat)
    })
    mfis_mat <- do.call(cbind, ss_clusterings)

    ### Extract negative population MFIs from UNSTAINED
    if (!is.null(sample_unstained)) {
        mfis_mat_unstained <- apply(sample_unstained, 2, median)

        # Make names easier to match
        # FS INT, FS TOF, SS INT, FL1 empty FITC, FL2 empty PE, FL3 empty ECD,
        # FL4 empty PC5.5, FL5 empty PC7, FL6 empty APC, FL7 empty AF700,
        # FL8 empty AA750, FL9 empty PB, FL10 empty KrO, TIME
        names(mfis_mat_unstained) <- sub("(FL[0-9]*).*", "\\1", names(mfis_mat_unstained))
        # FS INT, FS TOF, SS INT, FL1, FL2, FL3, FL4, FL5, FL6, FL7, FL8, FL9, FL10, TIME

        # Now replace the first row of mfis_mat (the minimum values) with the UNSTAINED values
        # select the columns of mfis_mat_unstained in the same order as mfis_mat
        # Then replace the first row of mfis_mat with the values of mfis_mat_unstained
        mfis_mat[1, ] <- mfis_mat_unstained[sub("(FL[0-9]*).*", "\\1", colnames(mfis_mat))]
    }
    return(mfis_mat)
}
