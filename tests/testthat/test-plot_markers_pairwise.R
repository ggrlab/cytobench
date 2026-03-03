test_that("Plot markers pairwise: Add mode-lines", {
    cn <- c(
        "FITC-A", "PE-A", "ECD-A"
    )
    tmpmat <- matrix(
        c(
            rnorm(500 * length(cn), mean = 1000, sd = 100),
            rnorm(1000 * length(cn), mean = 0, sd = 100)
        ),
        ncol = length(cn),
        byrow = TRUE
    ) %*% diag(c(3, 1, 1))
    tmpmat[, 1] <- tmpmat[, 1] + 1e5
    colnames(tmpmat) <- cn
    example_ff <- flowCore::flowFrame(exprs = tmpmat)
    tmpdir <- local_tempdir_time()

    #### Here I have to care that
    # 1) hdrcde is installed
    # 2) The mode lines are drawn in the "major" population (bottom-left at around x=0, y=0)

    p0 <- plot_markers_pairwise(
        example_ff,
        title_global = "Global title, e.g. sample ID",
        modelines = TRUE,
        kwargs_hdr = list(prob = 0.99)
    )
    pdf(file.path(tmpdir, "removeme.pdf"), height = 45, width = 45)
    print(p0)
    dev.off()

    pdf(file.path(tmpdir, "removeme2.pdf"), height = 5, width = 5)
    plot_markers_pairwise(
        example_ff,
        engine = "base", geom = "points",
        title_global = "Global title, e.g. sample ID",
        modelines = TRUE,
        kwargs_hdr = list(prob = 0.99)
    )
    dev.off()
    browser()
    testthat::expect_true(TRUE) # run through test
})
stop()

test_that("Plot markers pairwise performance", {
    cn <- c(
        "FITC-A", "PE-A", "ECD-A"
    )
    npoint <- 1e4
    example_ff <- flowCore::flowFrame(
        matrix(
            rnorm(npoint * length(cn), mean = 1000, sd = 100),
            ncol = length(cn),
            dimnames = list(NULL, cn)
        )
    )
    p0 <- plot_markers_pairwise(
        example_ff
    )
    tmpdir <- local_tempdir_time()
    times <- c()
    times <- c(times, Sys.time())
    pdf(file.path(tmpdir, "removeme.pdf"), height = 5, width = 5)
    print(p0)
    dev.off()
    times <- c(times, Sys.time())

    pdf(file.path(tmpdir, "removeme2.pdf"), height = 5, width = 5)
    plot_markers_pairwise(
        example_ff,
        engine = "base", geom = "points",
    )
    dev.off()
    times <- c(times, Sys.time())
    print(diff(times))
    testthat::expect_lt(times[1], times[2], label = "ggplot time", expected.label = "base time")
})

test_that("Plot markers pairwise performance_v2", {
    cn <- c(
        "FITC-A", "PE-A", "ECD-A"
    )
    cofactor_namedvec <- setNames(
        rep(5, length(cn)),
        cn
    )
    npoint <- 1e4
    example_ff <- flowCore::flowFrame(
        matrix(
            rnorm(npoint * length(cn), mean = 1000, sd = 100),
            ncol = length(cn),
            dimnames = list(NULL, cn)
        )
    )
    tmpdir <- local_tempdir_time()
    times <- c()
    times <- c(times, Sys.time())
    pdf(file.path(tmpdir, "removeme.pdf"), height = 5, width = 5)
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
    print(p0)
    dev.off()
    times <- c(times, Sys.time())

    pdf(file.path(tmpdir, "removeme.pdf"), height = 5, width = 5)
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
        geom = "hex"
    )
    print(p1)
    dev.off()
    times <- c(times, Sys.time())
    pdf(file.path(tmpdir, "removeme2.pdf"), height = 5, width = 5)
    plot_markers_pairwise(
        example_ff,
        cofactor_namedvec = cofactor_namedvec,
        special_cofactor_list = list(),
        transform_fun = function(x) {
            log10(x)
        },
        transform_fun_name = "log10",
        geom = "points",
        engine = "base",
        cex = 50,
    )
    dev.off()
    times <- c(times, Sys.time())
    pdf(file.path(tmpdir, "removeme2.pdf"), height = 5, width = 5)
    plot_markers_pairwise(
        example_ff,
        cofactor_namedvec = cofactor_namedvec,
        special_cofactor_list = list(),
        transform_fun = function(x) {
            log10(x)
        },
        transform_fun_name = "log10",
        geom = "points",
        engine = "base"
    )
    dev.off()
    times <- c(times, Sys.time())


    print(diff(times))
    testthat::expect_lt(
        object = times[1], label = "ggplot point time",
        expected = times[2], expected.label = "ggplot hex time"
    )
    testthat::expect_lt(
        object = times[2], label = "ggplot hex time",
        expected = times[3], expected.label = "base_cex5 time"
    )
    testthat::expect_lt(
        object = times[3], label = "base_cex25 time",
        expected = times[4], expected.label = "base time"
    )
})


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

    pdf(file.path(tmpdir, "removeme2.pdf"), height = 5, width = 5)
    plot_markers_pairwise(
        example_ff,
        engine = "base", geom = "points",
    )
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
        rep(5, length(cn)),
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


    pdf(file.path(tmpdir, "removeme2.pdf"), height = 5, width = 5)
    plot_markers_pairwise(
        example_ff,
        engine = "base", geom = "points",
        cofactor_namedvec = cofactor_namedvec,
        special_cofactor_list = list(),
        transform_fun = function(x) {
            log10(x)
        },
        transform_fun_name = "log10",
    )
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
        rep(5, length(cn)),
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
    times <- c()
    tmpdir <- local_tempdir_time()
    times <- c(times, Sys.time())
    pdf(file.path(tmpdir, "removeme.pdf"), height = 45, width = 45)
    print(p0)
    dev.off()
    times <- c(times, Sys.time())
    pdf(file.path(tmpdir, "removeme.pdf"), height = 45, width = 45)
    print(p1)
    dev.off()
    for (geom_x in c("points", "pointdensity")) {
        times <- c(times, Sys.time())
        pdf(file.path(tmpdir, "removeme.pdf"), height = 5, width = 5)
        plot_markers_pairwise(
            example_ff,
            cofactor_namedvec = cofactor_namedvec,
            special_cofactor_list = list(),
            transform_fun = function(x) {
                log10(x)
            },
            transform_fun_name = "log10",
            geom = geom_x,
            engine = "base",
            cex = 5,
        )
        dev.off()
    }
    times <- c(times, Sys.time())
    print(file.path(tmpdir, "removeme.pdf"))
    print(diff(times))


    testthat::expect_true(TRUE) # run through test
})


test_that("Plot markers pairwise title", {
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
        example_ff,
        title_global = "Global title, e.g. sample ID"
    )
    tmpdir <- local_tempdir_time()
    pdf(file.path(tmpdir, "removeme.pdf"), height = 45, width = 45)
    print(p0)
    dev.off()

    pdf(file.path(tmpdir, "removeme2.pdf"), height = 5, width = 5)
    plot_markers_pairwise(
        example_ff,
        engine = "base", geom = "points",
        title_global = "Global title, e.g. sample ID"
    )
    dev.off()

    testthat::expect_true(TRUE) # run through test
})
devtools::load_all()
test_that("Plot markers pairwise with density", {
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
    tmpdir <- local_tempdir_time()

    pdf(file.path(tmpdir, "removeme2.pdf"), height = 5, width = 5)
    plot_markers_pairwise(
        example_ff,
        engine = "base", geom = "points",
        title_global = "Global title, e.g. sample ID",
        diagonal_densityplot = TRUE
    )
    dev.off()

    testthat::expect_true(TRUE) # run through test
})


test_that("Check layoutcreation", {
    cn <- c("FITC-A", "PE-A", "ECD-A")
    expected_nodiag <- matrix(
        c(0, 1, 2, 3, 4, 5, 6, 0, 7),
        nrow = 3,
        dimnames = list(
            c("namecol", "PE-A", "ECD-A"),
            c("namecol", "FITC-A", "PE-A")
        )
    )
    expected_diag <- matrix(
        c(0:3, 4:7, 8, 0, 9, 10, 11, 0, 0, 12),
        nrow = 4,
        dimnames = list(
            c("namecol", "FITC-A", "PE-A", "ECD-A"),
            c("namecol", "FITC-A", "PE-A", "ECD-A")
        )
    )
    layout_nodiag <- create_layout(cn)
    testthat::expect_equal(layout_nodiag, expected = expected_nodiag)

    layout_diag <- create_layout(cn, diagonal_densityplot = TRUE)
    testthat::expect_equal(layout_diag, expected = expected_diag)
})

testthat::test_that("Expect warning with base engine and NO geom", {
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
    tmpdir <- local_tempdir_time()
    testthat::expect_warning(
        plot_markers_pairwise(
            example_ff,
            engine = "base",
            title_global = "Global title, e.g. sample ID"
        ),
        "geom = 'hex' not implemented in base engine; using 'points' instead"
    )
})
