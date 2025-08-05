test_that("Export files as csv", {
    dt_simulated <- simulate_fs(5, flowcore = FALSE)

    Sys.sleep(.01) # Ensure different timestamps
    exported_files <- export_csv(dt_simulated, outdir = local_tempdir_time())
    testthat::expect_equal(names(exported_files), names(dt_simulated))

    read_exported_files <- lapply(exported_files, data.table::fread)
    for (i in seq_len(length(read_exported_files))) {
        testthat::expect_equal(
            read_exported_files[[i]],
            dt_simulated[[i]]
        )
    }
})

test_that("Export files as feather", {
    dt_simulated <- simulate_fs(5, flowcore = FALSE)

    Sys.sleep(.01) # Ensure different timestamps
    exported_files <- export_feather(dt_simulated, outdir = local_tempdir_time())
    testthat::expect_equal(names(exported_files), names(dt_simulated))


    read_exported_files <- lapply(exported_files, function(x) {
        data.table::as.data.table(feather::read_feather(x))
    })
    for (i in seq_len(length(read_exported_files))) {
        testthat::expect_equal(
            read_exported_files[[i]],
            dt_simulated[[i]]
        )
    }
})
