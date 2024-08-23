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
    return(fSOM)
}
test_that("FlowSOM outliers", {
    fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
    ff <- flowCore::read.FCS(fileName)
    fsom <- do_flowsom()
    outliers <- flowSOM_is.outlier(fsom)
    testthat::expect_equal(sum(outliers), 209)
    testthat::expect_length(outliers, flowCore::nrow(ff))
})



test_that("FlowSOM outliers", {
    ff_example <- example_processed()
    fsom <- do_flowsom(ff_example)
    res <- flowSOM_predict(flowsom_result = fsom, flowset = flowCore::flowSet(ff_example))

    outliers <- flowSOM_is.outlier(fsom)
    testthat::expect_equal(outliers, res$cell_outlier)
    testthat::expect_equal(sum(outliers), 209)
    testthat::expect_equal(sum(res$cell_outlier), 209)
})
