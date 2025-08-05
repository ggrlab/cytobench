test_that("Transform data", {
    cn <- c(
        "FS INT", "FS TOF", "SS INT",
        "FL1 CD45RA FITC", "FL2 CCR7 PE", "FL3 CD28 ECD", "FL4 PD1 PC5.5",
        "FL5 CD27 PC7", "FL6 CD4 APC", "FL7 CD8 AF700", "FL8 CD3 AA750",
        "FL9 CD57 PB", "FL10 CD45 KrO", "TIME"
    )
    dt_simulated <- simulate_fs(3, flowcore = FALSE, columns = cn)
    dt_simulated_backup <- lapply(dt_simulated, data.table::data.table)

    transformed_no_inplace <- transform_dts(dt_simulated, inplace = FALSE)
    # transformed must be different from csv_read
    testthat::expect_true(
        all((transformed_no_inplace[[1]] != dt_simulated[[1]])[, -c(1:3, 14)])
    )
    # dt_simulated must be the same as dt_simulated_backup
    testthat::expect_true(all((dt_simulated[[1]] == dt_simulated_backup[[1]])))


    transformed_inplace <- transform_dts(dt_simulated, inplace = TRUE)
    # If inplace was done, transformed_inplace must be the same as dt_simulated
    testthat::expect_true(
        all((transformed_no_inplace[[1]] == dt_simulated[[1]]))
    )
    # and dt_simulated must be different from dt_simulated_backup

    # dt_simulated must be the same as dt_simulated_backup
    testthat::expect_true(
        all((dt_simulated[[1]] != dt_simulated_backup[[1]])[, -c(1:3, 14)])
    )
})


test_that("Transform data different function", {
    cn <- c(
        "FS INT", "FS TOF", "SS INT",
        "FL1 CD45RA FITC", "FL2 CCR7 PE", "FL3 CD28 ECD", "FL4 PD1 PC5.5",
        "FL5 CD27 PC7", "FL6 CD4 APC", "FL7 CD8 AF700", "FL8 CD3 AA750",
        "FL9 CD57 PB", "FL10 CD45 KrO", "TIME"
    )
    dt_simulated <- simulate_fs(3, flowcore = FALSE, columns = cn)

    transformed_weird_logarithm_do_not_use <- transform_dts(
        dt_simulated,
        transform_fun = function(x_values, scaling_factor) {
            log(abs(x_values) / scaling_factor)
        },
        inplace = FALSE
    )
    transformed_weird_logarithm_do_not_use <- transform_dts(
        dt_simulated,
        transform_fun = function(x_values, scaling_factor) {
            asinh(x_values / scaling_factor)
        },
        inplace = FALSE
    )
    testthat::expect_true(TRUE) # Run through test
})
