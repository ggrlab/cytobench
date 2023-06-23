devtools::load_all()
library(testthat)
test_that("Export files as csv", {
    # Read in the example data
    csv_paths <- list.files(
        system.file("extdata", "raw_csv", package = "cytobench"),
        full.names = TRUE, recursive = TRUE,
        pattern = "*.csv$"
    )
    csv_read <- lapply(csv_paths, data.table::fread)
    names(csv_read) <- sub(".*raw_csv/", "", csv_paths)

    exported_files <- export_csv(csv_read, outdir = tempdir())
    testthat::expect_equal(names(exported_files), names(csv_read))

    read_exported_files <- lapply(exported_files, data.table::fread)
    for (i in seq_len(length(read_exported_files))) {
        testthat::expect_equal(
            read_exported_files[[i]],
            csv_read[[i]]
        )
    }
})

test_that("Export files as feather", {
    # Read in the example data
    csv_paths <- list.files(
        system.file("extdata", "raw_csv", package = "cytobench"),
        full.names = TRUE, recursive = TRUE,
        pattern = "*.csv$"
    )
    csv_read <- lapply(csv_paths, data.table::fread)
    names(csv_read) <- sub(".*raw_csv/", "", csv_paths)

    exported_files <- export_feather(csv_read, outdir = tempdir())
    testthat::expect_equal(names(exported_files), names(csv_read))


    read_exported_files <- lapply(exported_files, function(x) {
        data.table::as.data.table(feather::read_feather(x))
    })
    for (i in seq_len(length(read_exported_files))) {
        testthat::expect_equal(
            read_exported_files[[i]],
            csv_read[[i]]
        )
    }
})
