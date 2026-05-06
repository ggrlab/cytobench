devtools::load_all()

test_that("Test gating cells and extraction", {
    gs <- flowWorkspace::load_gs(testthat::test_path("testdata", "HKP P0Y V01 Basic 001"))
    flowWorkspace::gs_pop_add(gs, flowCore::rectangleGate("FSC_A" = c(0, 10)))
    flowWorkspace::recompute(gs)


    a1 <- extract_gatingset_cells(gs, gatenames = "/Singlets/CD45+")
    testthat::expect_true(inherits(a1[["extracted"]][[1]], "flowFrame"))
    a2 <- extract_gatingset_cells(gs, gatenames = "/Singlets/CD45+", return_x = "cytoset")
    testthat::expect_true(inherits(a2[["extracted"]][[1]], "cytoframe"))
    a3 <- extract_gatingset_cells(gs, gatenames = "/Singlets/CD45+", return_x = "gatinghierarchy")
    testthat::expect_true(inherits(a3[["extracted"]], "GatingSet"))
})


test_that("Test gating cells and returning a data.table with the gate ids!", {
    gs <- flowWorkspace::load_gs(testthat::test_path("testdata", "HKP P0Y V01 Basic 001"))
    flowWorkspace::gs_pop_add(gs, flowCore::rectangleGate("FSC_A" = c(0, 10)))
    flowWorkspace::recompute(gs)

    assigned_cells <- assign_gatingset_cells(gs)
    testthat::expect_equal(
        colnames(assigned_cells),
        flowWorkspace::gs_get_pop_paths(gs),
        info = "Gate paths should be the column names of the assigned cells data.table. In this particular case, no cells are present, so no rows are expected."
    )

    assigned_cells_withexprs <- assign_gatingset_cells(gs, return_exprs = TRUE)
    testthat::expect_equal(
        colnames(assigned_cells_withexprs),
        c(colnames(flowWorkspace::gs_get_cytoframe(gs, 1)), flowWorkspace::gs_get_pop_paths(gs)),
        info = "Column names should include both gate paths and expression column names."
    )


    assigned_cells_specific <- assign_gatingset_cells(gs, gatenames = c("/Singlets/CD45+"))
    testthat::expect_equal(
        colnames(assigned_cells_specific),
        "/Singlets/CD45+",
        info = "Column names should match the specified gatenames, resolving numeric indices to gate paths."
    )


    assigned_cells_specificnumbered <- assign_gatingset_cells(gs, gatenames = c(5))
    testthat::expect_equal(
        colnames(assigned_cells_specificnumbered),
        flowWorkspace::gs_get_pop_paths(gs)[5],
        info = "Column names should match the specified gatenames, resolving numeric indices to gate paths."
    )
    testthat::expect_equal(
        colnames(assigned_cells_specificnumbered),
        "/Singlets/CD45+/D/H",
        info = "The 5th gate path should be '/Singlets/CD45+/D/H'."
    )
})
