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
        engine = "base"
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
        engine = "base"
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
        engine = "base",
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
