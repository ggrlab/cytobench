# devtools::load_all()

test_that("Plot markers pairwise", {
    example_ff <- flowCore::flowFrame(
        matrix(
            rnorm(1000 * 10, mean = 1000, sd = 100),
            ncol = 10,
            dimnames = list(NULL, c(
                "FITC-A", "PE-A", "ECD-A", "PC5.5-A", "PC7-A",
                "APC-A", "AF700-A", "AA750-A", "PB-A", "KrO-A"
            ))
        )
    )
    p0 <- plot_markers_pairwise(
        example_ff
    )
    pdf("removeme.pdf", height = 45, width = 45)
    print(p0)
    dev.off()
})


test_that("Plot markers pairwise: Cofactors", {
    cn <- c(
        "FITC-A", "PE-A", "ECD-A", "PC5.5-A", "PC7-A",
        "APC-A", "AF700-A", "AA750-A", "PB-A", "KrO-A"
    )
    example_ff <- flowCore::flowFrame(
        matrix(
            rnorm(1000 * 10, mean = 1000, sd = 100),
            ncol = 10,
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
        transform_fun = function(x){
            log10(x)
        },
        transform_fun_name = "log10",
    )
    pdf("removeme.pdf", height = 45, width = 45)
    print(p0)
    dev.off()
})

test_that("Plot markers pairwise: Points", {
    cn <- c(
        "FITC-A", "PE-A", "ECD-A", "PC5.5-A", "PC7-A",
        "APC-A", "AF700-A", "AA750-A", "PB-A", "KrO-A"
    )
    example_ff <- flowCore::flowFrame(
        matrix(
            rnorm(1000 * 10, mean = 1000, sd = 100),
            ncol = 10,
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
        transform_fun = function(x){
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
        transform_fun = function(x){
            log10(x)
        },
        transform_fun_name = "log10",
        geom = "pointdensity",
        bins = 250
    )
    # browser()
    # pdf("removeme.pdf", height = 45, width = 45)
    # print(p0)
    # print(p1)
    # dev.off()
})
