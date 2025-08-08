test_that("Plot markers pairwise", {
    cn <- c(
        "FITC-A", "PE-A", "ECD-A"
    )
    example_ff <- flowCore::flowFrame(
        matrix(
            rnorm(1000 * length(cn), mean = 1000, sd = 100),
            ncol = length(cn),
            dimnames = list(NULL, cn)
        )
    )
    p0 <- plot_markers_pairwise(
        example_ff
    )
    tmpdir <- local_tempdir_time()
    pdf(file.path(tmpdir, "removeme.pdf"), height = 45, width = 45)
    print(p0)
    dev.off()
    testthat::expect_true(TRUE) # run through test
})


test_that("Plot markers pairwise: Cofactors", {
    cn <- c(
        "FITC-A", "PE-A", "ECD-A"
    )
    example_ff <- flowCore::flowFrame(
        matrix(
            rnorm(1000 * length(cn), mean = 1000, sd = 100),
            ncol = length(cn),
            dimnames = list(NULL, cn)
        )
    )
    cofactor_namedvec <- setNames(
        c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5),
        cn
    )
    p0 <- plot_markers_pairwise(
        example_ff,
        cofactor_namedvec = cofactor_namedvec,
        special_cofactor_list = list(
            "FITC-A_PE-A" = c(1000, 1000)
        ),
        transform_fun = function(x) {
            log10(x)
        },
        transform_fun_name = "log10",
    )
    tmpdir <- local_tempdir_time()
    pdf(file.path(tmpdir, "removeme.pdf"), height = 45, width = 45)
    print(p0)
    dev.off()
    testthat::expect_true(TRUE) # run through test
})

test_that("Plot markers pairwise: Points", {
    cn <- c(
        "FITC-A", "PE-A", "ECD-A"
    )
    example_ff <- flowCore::flowFrame(
        matrix(
            rnorm(1000 * length(cn), mean = 1000, sd = 100),
            ncol = length(cn),
            dimnames = list(NULL, cn)
        )
    )
    cofactor_namedvec <- setNames(
        c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5),
        cn
    )
    p0 <- plot_markers_pairwise(
        example_ff,
        cofactor_namedvec = cofactor_namedvec,
        special_cofactor_list = list(
            "FITC-A_PE-A" = c(1000, 1000)
        ),
        transform_fun = function(x) {
            log10(x)
        },
        transform_fun_name = "log10",
        geom = "points"
    )
    p1 <- plot_markers_pairwise(
        example_ff,
        cofactor_namedvec = cofactor_namedvec,
        special_cofactor_list = list(
            "FITC-A_PE-A" = c(1000, 1000)
        ),
        transform_fun = function(x) {
            log10(x)
        },
        transform_fun_name = "log10",
        geom = "pointdensity",
        bins = 250
    )
    tmpdir <- local_tempdir_time()
    pdf(file.path(tmpdir, "removeme.pdf"), height = 45, width = 45)
    print(p0)
    print(p1)
    dev.off()
    testthat::expect_true(TRUE) # run through test
})
