# devtools::load_all()

test_that("FlowSOM outliers", {
    fileName <- system.file("extdata", "68983.fcs", package = "FlowSOM")
    ff <- flowCore::read.FCS(fileName)
    # Compensation
    comp <- flowCore::keyword(ff)[["SPILL"]]
    ff <- flowWorkspace::compensate(ff, comp)
    # Transformation
    transformList <- flowCore::estimateLogicle(ff, channels = colnames(comp))
    ff <- flowWorkspace::transform(ff, transformList)

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

    outliers <- flowSOM_is.outlier(fSOM)
    testthat::expect_equal(sum(outliers), 209)
    testthat::expect_length(outliers, flowCore::nrow(ff))
})
