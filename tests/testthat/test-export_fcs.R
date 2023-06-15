test_that("Export fcs files", {
    # Read in the example data
    csv_paths <- list.files(
        system.file("extdata", "raw_csv", package = "cytobench"),
        full.names = TRUE, recursive = TRUE, 
        pattern="*.csv$"
    )
    csv_read <- lapply(csv_paths, data.table::fread)
    names(csv_read) <- csv_paths

    tmpdir <- tempdir()
    # Write the fcs files
    export_fcs(
        csv_read,
        outdir = file.path(tmpdir, "fcs_noshift"),
        safety_scaling = 1.25, 
        safety_shift = 0, 
        use.names = FALSE, 
        new_colnames_to = "description", 
        new_colnames = 1
    )
})
