devtools::load_all()
test_that("Transform data", {
    # Read in the example data
    csv_paths <- list.files(
        system.file("extdata", "raw_csv", package = "cytobench"),
        full.names = TRUE, recursive = TRUE,
        pattern = "*.csv$"
    )
    csv_read <- lapply(csv_paths, data.table::fread)
    csv_read_backup <- lapply(csv_paths, data.table::fread)
    names(csv_read) <- sub(".*raw_csv/", "", csv_paths)
    names(csv_read_backup) <- sub(".*raw_csv/", "", csv_paths)

    transformed_no_inplace <- transform(csv_read[1:3], inplace = FALSE)
    # transformed must be different from csv_read
    testthat::expect_true(
        all((transformed_no_inplace[[1]] != csv_read[[1]])[, -c(1:3, 14)])
    )
    # csv_read must be the same as csv_read_backup
    testthat::expect_true(all((csv_read[[1]] == csv_read_backup[[1]])))


    transformed_inplace <- transform(csv_read[1:3], inplace = TRUE)
    # If inplace was done, transformed_inplace must be the same as csv_read
    testthat::expect_true(
        all((transformed_no_inplace[[1]] == csv_read[[1]]))
    )
    # and csv_read must be different from csv_read_backup

    # csv_read must be the same as csv_read_backup
    testthat::expect_true(
        all((csv_read[[1]] != csv_read_backup[[1]])[, -c(1:3, 14)])
    )
})
test_that("Transform data different function", {
    # Read in the example data
    csv_paths <- list.files(
        system.file("extdata", "raw_csv", package = "cytobench"),
        full.names = TRUE, recursive = TRUE,
        pattern = "*.csv$"
    )
    csv_read <- lapply(csv_paths, data.table::fread)
    names(csv_read) <- sub(".*raw_csv/", "", csv_paths)

    transformed_weird_logarithm_do_not_use <- transform(
        csv_read[1:3],
        transform_fun = function(x_values, scaling_factor) {
            log(abs(x_values) / scaling_factor)
        },
        inplace = FALSE
    )
    transformed_weird_logarithm_do_not_use <- transform(
        csv_read[1:3],
        transform_fun = function(x_values, scaling_factor) {
            asinh(x_values / scaling_factor)
        },
        inplace = FALSE
    )
})
