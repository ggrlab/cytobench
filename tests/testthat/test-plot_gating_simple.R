test_that("Test plot gating simple", {
    
    gs <- flowWorkspace::load_gs(testthat::test_path("testdata", "HKP P0Y V01 Basic 001"))
    flowWorkspace::gs_pop_add(gs, flowCore::rectangleGate("FSC_A" = c(0, 10)))
    flowWorkspace::recompute(gs)
    plots <- plot_gating_simple(gs, facet = FALSE)
    pdf("removeme.pdf")
    print(plots)
    dev.off()
    testthat::expect_true(TRUE)
})
