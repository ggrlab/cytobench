test_that("Write an FCS file", {
    ff <- simulate_ff()
    tmp <- write_memory_FCS(ff)
    reread_ff <- flowCore::read.FCS(tmp)
    testthat::expect_false(identical(reread_ff, ff)) # they SHOULD be different objects!
    testthat::expect_equal(nrow(reread_ff), nrow(ff))
    testthat::expect_true(max(abs(flowCore::exprs(ff) - flowCore::exprs(reread_ff))) < 1e-8)
})
