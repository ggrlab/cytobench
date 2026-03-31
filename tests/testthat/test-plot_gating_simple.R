test_that("Test plot gating simple", {
    gs <- flowWorkspace::load_gs(testthat::test_path("testdata", "HKP P0Y V01 Basic 001"))
    flowWorkspace::gs_pop_add(gs, flowCore::rectangleGate("FSC_A" = c(0, 10)))
    flowWorkspace::recompute(gs)
    plots <- plot_gating_simple(gs, facet = FALSE)
    # Expect warning: "this is a simple plot, only showing the first 2 levels of the gating hierarchy"

    pdf("removeme.pdf")
    w <- testthat::capture_warnings(print(plots))
    dev.off()
    #  No shared levels found between `names(values)` of the manual scale and the data's colour values.
    testthat::expect_true(all(grepl("No shared levels found between .* of the manual scale and the data's .* values.", w, fixed = FALSE)))
    # The warning occurs because there are NO cells in gs, so the plot is empty and ggplot complains
    testthat::expect_true(TRUE)
})
