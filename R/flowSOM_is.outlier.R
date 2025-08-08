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
#' @keywords flowsom
#' @examples
#' ff_example <- example_processed()
#' fsom <- do_flowsom_TESTING(ff_example)
#' outliers <- flowSOM_is.outlier(fsom)
flowSOM_is.outlier <- function(fsom,
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

    # The following is a copy of the FlowSOM::TestOutliers() function,
    # the comments were mostly added by me.
    fsom <- FlowSOM::UpdateFlowSOM(fsom)
    # If no reference FlowSOM object is provided, use the original fsom object
    if (is.null(fsomReference)) {
        fsomReference <- fsom
    }

    # Get cluster assignments for the reference FlowSOM object as factor
    referenceClusters <- factor(FlowSOM::GetClusters(fsomReference),
        levels = seq_len(fsomReference$map$nNodes)
    )

    # Get cluster assignments for the current FlowSOM object as factor
    clusters <- factor(FlowSOM::GetClusters(fsom),
        levels = seq_len(fsom$map$nNodes)
    )

    # Distance in high-dimensional space [is written in fsomReference$map$mapping[, 2]]
    # Calculate the median values within each reference cluster
    medians <- tapply(
        fsomReference$map$mapping[, 2],
        referenceClusters,
        stats::median
    )

    # Calculate the median absolute deviations (MAD) within each reference cluster
    mads <- tapply(
        fsomReference$map$mapping[, 2],
        referenceClusters,
        stats::mad
    )

    # Determine the threshold within outliers based on the allowed deviation in terms of MAD
    thresholds <- medians + madAllowed * mads

    # Identify outliers by comparing the mapping values to the calculated thresholds
    # of the cell's assigned cluster
    outliers <- (fsom$map$mapping[, 2] > thresholds[clusters])
    return(unname(outliers))
}
