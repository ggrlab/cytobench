test_that("Export fcs files, new colnames into description", {
    # Create the example data
    original_cn <- c(
        "FITC-A",
        "PE-A",
        "ECD-A",
        "PC5.5-A",
        "PC7-A",
        "APC-A",
        "AF700-A",
        "AA750-A",
        "PB-A",
        "KrO-A"
    )
    dt_simulated <- simulate_fs(5, flowcore = FALSE, columns = original_cn)
    new_cn <- c(
        "FL1 CD45RA FITC", "FL2 CCR7 PE",
        "FL3 CD28 ECD", "FL4 PD1 PC5.5", "FL5 CD27 PC7", "FL6 CD4 APC",
        "FL7 CD8 AF700", "FL8 CD3 AA750", "FL9 CD57 PB", "FL10 CD45 KrO"
    )
    dt_simulated_newcolnames <- simulate_fs(5, flowcore = FALSE, columns = new_cn)
    # I actively test if spaces work here.
    names(dt_simulated_newcolnames) <- paste0("newcol - ", names(dt_simulated_newcolnames))
    both <- c(
        dt_simulated,
        dt_simulated_newcolnames
    )

    tmpdir <- local_tempdir_time()
    # Write the fcs files
    extreme_fcs <- export_fcs(
        both,
        outdir = file.path(tmpdir, "fcs_noshift"),
        safety_scaling = 1.25,
        safety_shift = 0,
        use.names = FALSE,
        new_colnames_to = "description",
        new_colnames = "simsample_1"
    )

    # validate the fcs
    validate <- function(expected_cn = original_cn) {
        ## read all fcs
        fcs_paths <- list.files(
            file.path(tmpdir, "fcs_noshift"),
            full.names = TRUE, recursive = TRUE,
            pattern = "*.fcs$"
        )
        fcs_read <- lapply(fcs_paths, flowCore::read.FCS)
        marker_names <- lapply(fcs_read, flowCore::markernames)
        for (sample_x_colnames in marker_names) {
            expect_identical(
                names(sample_x_colnames),
                expected_cn
            )
        }
    }
    validate()
    unlink(file.path(tmpdir, "fcs_noshift"), recursive = TRUE)
    extreme_fcs_both <- export_fcs(
        both,
        outdir = file.path(tmpdir, "fcs_noshift"),
        safety_scaling = 1.25,
        safety_shift = 0,
        use.names = FALSE,
        new_colnames_to = "description",
        new_colnames = "newcol - simsample_5"
    )
    validate(expected_cn = new_cn)
    extreme_fcs_both_v2 <- export_fcs(
        both,
        outdir = file.path(tmpdir, "fcs_noshift"),
        safety_scaling = 1.25,
        safety_shift = 0,
        use.names = FALSE,
        new_colnames_to = "description"
    )
    validate(expected_cn = original_cn)

    extreme_fcs_both_v3 <- export_fcs(
        both,
        outdir = file.path(tmpdir, "fcs_noshift"),
        safety_scaling = 1.25,
        safety_shift = 0,
        use.names = FALSE,
        new_colnames_to = "description",
        new_colnames = 1
    )
    validate(expected_cn = original_cn)
})

test_that("Export fcs files, shifts", {
    dt_simulated <- simulate_fs(3, flowcore = FALSE)

    tmpdir <- local_tempdir_time()
    # Write the fcs files
    extreme_fcs <- export_fcs(
        dt_simulated,
        outdir = file.path(tmpdir, "fcs_noshift_shift_0"),
        safety_scaling = 1.25,
        safety_shift = 0,
        use.names = FALSE,
        new_colnames = 1
    )
    extreme_fcs <- export_fcs(
        dt_simulated,
        outdir = file.path(tmpdir, "fcs_noshift_shift_10"),
        safety_scaling = 10.25,
        safety_shift = 10,
        use.names = FALSE,
        new_colnames = 1
    )
    extreme_fcs <- export_fcs(
        dt_simulated,
        outdir = file.path(tmpdir, "fcs_noshift_shift_neg10"),
        safety_scaling = 0.25,
        safety_shift = -10,
        use.names = FALSE,
        new_colnames = 1
    )

    # validate the fcs
    ## read all fcs
    read_fcs <- lapply(
        list(
            "0" = file.path(tmpdir, "fcs_noshift_shift_0"),
            "10" = file.path(tmpdir, "fcs_noshift_shift_10"),
            "-10" = file.path(tmpdir, "fcs_noshift_shift_neg10")
        ),
        function(x) {
            fcs_paths <- list.files(
                x,
                full.names = TRUE, recursive = TRUE,
                pattern = "*.fcs$"
            )
            fcs_read <- lapply(fcs_paths, flowCore::read.FCS)
        }
    )
    for (i in 1:length(read_fcs[[1]])) {
        d_0 <- flowCore::exprs(read_fcs[["0"]][[i]])
        d_10 <- flowCore::exprs(read_fcs[["10"]][[i]])
        d_neg10 <- flowCore::exprs(read_fcs[["-10"]][[i]])

        # remove attributes
        attr(d_0, "ranges") <- NULL
        attr(d_10, "ranges") <- NULL
        attr(d_neg10, "ranges") <- NULL

        expect_equal(
            d_0,
            d_10 - 10,
            tolerance = 1e-5
        )
        expect_identical(
            d_0,
            d_neg10 + 10,
            tolerance = 1e-5
        )
    }
})

test_that("Export fcs files, existing template", {
    dt_simulated <- simulate_fs(3, flowcore = FALSE)

    tmpdir <- local_tempdir_time()
    # Write the fcs files
    extreme_fcs <- export_fcs(
        dt_simulated,
        outdir = file.path(tmpdir, "extreme_fcs_original"),
        safety_scaling = 1.25,
        safety_shift = 0,
        use.names = FALSE,
        new_colnames = 1
    )
    extreme_fcs_v2 <- export_fcs(
        dt_simulated,
        outdir = file.path(tmpdir, "extreme_fcs_template"),
        safety_scaling = 10.25,
        safety_shift = 0,
        extreme_template = extreme_fcs
    )

    # validate the fcs
    ## read all fcs
    read_fcs <- lapply(
        list(
            "original" = file.path(tmpdir, "extreme_fcs_original"),
            "template" = file.path(tmpdir, "extreme_fcs_template")
        ),
        function(x) {
            fcs_paths <- list.files(
                x,
                full.names = TRUE, recursive = TRUE,
                pattern = "*.fcs$"
            )
            fcs_read <- lapply(fcs_paths, flowCore::read.FCS)
        }
    )
    for (i in 1:length(read_fcs[[1]])) {
        expect_equal(
            flowCore::exprs(read_fcs[["original"]][[i]]),
            flowCore::exprs(read_fcs[["template"]][[i]])
        )
    }
})
