test_that("Plotting densities", {
    fs <- simulate_fs(n_samples = 3, ncells = 250, columns = c("CD4", "CD8"))
    fs_dt <- flowCore::fsApply(fs, function(x) data.table::as.data.table(flowCore::exprs(x))) |>
        data.table::rbindlist(idcol = "File")
    testthat::expect_error(
        plot_densities(fs_dt),
        regexp = "Expression values after melting cannot be character. Did you set groupings and markers correctly?"
    )
    plot_densities(fs_dt, markers = c("CD4", "CD8"))
    plot_densities(fs_dt, column_color = "File")

    # Explicitely calculating densities first, then plotting
    densities <- calc_densities(
        dt = fs_dt,
        groupings = "File"
    )
    plot_densities(densities, column_color = "File")
    plot_densities(densities)

    testthat::expect_error(
        densities <- calc_densities(dt = fs_dt, groupings = NULL),
        regexp = "Expression values after melting cannot be character. Did you set groupings and markers correctly?"
    )
    densities <- calc_densities(
        dt = fs_dt,
        groupings = NULL,
        markers = c("CD4", "CD8")
    )
    plot_densities(densities)
})

test_that("Plotting densities BASEPLOT", {
    fs <- simulate_fs(n_samples = 3, ncells = 250, columns = c("CD4", "CD8"))
    fs_dt <- flowCore::fsApply(fs, function(x) data.table::as.data.table(flowCore::exprs(x))) |>
        data.table::rbindlist(idcol = "File")

    plot_densities_base(fs_dt, markers = c("CD4"))
    plot_densities_base(fs_dt, markers = c("CD4", "CD8"))
    plot_densities_base(fs_dt, markers = c("CD4", "CD4"))
    plot_densities_base(fs_dt, markers = c("CD4"), column_color = "File")
    testthat::expect_true(TRUE)
})
