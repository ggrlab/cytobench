#' @title Is a cell an outlier?
#' @description
#' Test if any cells are too far from their cluster centers.
#' This code is based on the FlowSOM package, FlowSOM::TestOutliers().
#' Here I do the exact same thing as in the FlowSOM package, but I return if the cell is an outlier or not.
#'
#' In contrast, the FlowSOM package returns summaries of that, OR a per-channel outlier assignment.
#' This channelwise outlier assignment is different from the per-default reported outliers in the package itself.
#'
#' @param fsom  FlowSOM object
#' @param madAllowed Number of median absolute deviations allowed. Default = 4.
#' @param fsomReference FlowSOM object to use as reference. If NULL (default),
#'                       the original fsom object is used.
#' @param channels If channels are given, return the result of FlowSOM::TestOutliers()
#' @export
flowSOM_outliers <- function(fsom,
                             madAllowed = 4,
                             fsomReference = NULL,
                             channels = NULL) {
    # https://github.com/SofieVG/FlowSOM/blob/master/R/3_buildMST.R#L239
    if (!all(is.null(channels))) {
        return(
            FlowSOM::TestOutliers(
                fsom = fsom,
                madAllowed = madAllowed,
                fsomReference = fsomReference,
                channels = channels
            )
        )
    }

    # The following is a copy of the FlowSOM::TestOutliers() function

    fsom <- FlowSOM::UpdateFlowSOM(fsom)
    if (is.null(fsomReference)) {
        fsomReference <- fsom
    } else {
        fsomReference <- FlowSOM::UpdateFlowSOM(fsomReference)
    }

    referenceClusters <- factor(FlowSOM::GetClusters(fsomReference),
        levels = seq_len(fsomReference$map$nNodes)
    )
    clusters <- factor(FlowSOM::GetClusters(fsom),
        levels = seq_len(fsom$map$nNodes)
    )

    # Distance in high-dimensional space
    medians <- tapply(
        fsomReference$map$mapping[, 2],
        referenceClusters,
        stats::median
    )

    mads <- tapply(
        fsomReference$map$mapping[, 2],
        referenceClusters,
        stats::mad
    )

    thresholds <- medians + madAllowed * mads
    browser()
    outliers <- (fsom$map$mapping[, 2] > thresholds[clusters])
}
