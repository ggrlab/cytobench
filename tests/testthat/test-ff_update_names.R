devtools::load_all()
mk_test_ff <- function() {
    exprs <- matrix(c(1, 2, 3, 4), ncol = 2)
    colnames(exprs) <- c("A", "B")
    ff <- flowCore::flowFrame(exprs)
    flowCore::markernames(ff) <- structure(c("oldA", "oldB"),
        names = c("A", "B")
    )
    ff
}
mk_test_ff_spill_named <- function() {
    ff <- mk_test_ff()

    spillmat <- matrix(c(1, 0.1, 0.1, 1), ncol = 2)
    colnames(spillmat) <- c("A", "B")
    flowCore::keyword(ff)[["$SPILLOVER"]] <- spillmat
    flowCore::keyword(ff)[["$P1N"]] <- "A"
    flowCore::keyword(ff)[["$P2N"]] <- "B"

    ff
}

mk_test_markermap <- function(colnames_new = c("X", "Y")) {
    data.frame(
        name = c("A", "B"),
        colnames = colnames_new,
        description_channel.flurochrome.marker.awh = c("CD3", "CD4"),
        stringsAsFactors = FALSE
    )
}


test_that("ff_update_names renames channels and marker names", {
    ff <- mk_test_ff_spill_named()
    map <- mk_test_markermap()

    spillovermat <- flowCore::spillover(ff)
    out <- cytobench:::ff_update_names(flowframe = ff, markermap = map)

    expect_identical(flowCore::colnames(out), c("X", "Y"))
    expect_identical(unname(flowCore::markernames(out)), c("CD3", "CD4"))
    expect_identical(names(flowCore::markernames(out)), c("X", "Y"))
})

test_that("ff_update_names renames channels and marker names", {
    ff <- mk_test_ff()
    map <- mk_test_markermap()

    testthat::expect_warning(
        out <- cytobench:::ff_update_names(flowframe = ff, markermap = map),
        "Could not extract spillover matrices: No spillover matrix stored in that sample"
    )

    expect_identical(flowCore::colnames(out), c("X", "Y"))
    expect_identical(unname(flowCore::markernames(out)), c("CD3", "CD4"))
    expect_identical(names(flowCore::markernames(out)), c("X", "Y"))
})

test_that("ff_update_names errors on missing map columns", {
    ff <- mk_test_ff()
    bad_map <- data.frame(name = c("A", "B"), colnames = c("X", "Y"))

    expect_error(
        cytobench:::ff_update_names(flowframe = ff, markermap = bad_map),
        "missing required columns"
    )
})

test_that("ff_update_names errors when channels are not mappable", {
    ff <- mk_test_ff()
    map <- mk_test_markermap()
    map$name <- c("A", "C")

    expect_error(
        cytobench:::ff_update_names(flowframe = ff, markermap = map),
        "missing from `markermap\\$name`"
    )
})

test_that("ff_update_names errors on duplicated target channel names", {
    ff <- mk_test_ff()
    map <- mk_test_markermap(colnames_new = c("X", "X"))

    expect_error(
        cytobench:::ff_update_names(flowframe = ff, markermap = map),
        "contain NA or duplicates"
    )
})

test_that("get_new_parameternames maps indexed spillover columns via $PnN", {
    ff <- mk_test_ff()
    flowCore::keyword(ff)[["$P1N"]] <- "A"
    flowCore::keyword(ff)[["$P2N"]] <- "B"
    spillmat <- matrix(c(1, 0.1, 0.1, 1), ncol = 2)
    colnames(spillmat) <- c("1", "2")

    out <- cytobench:::get_new_parameternames(
        ff = ff,
        spillmat = spillmat,
        new_channel_names = c("X", "Y"),
        markermap_names = c("A", "B")
    )

    expect_identical(out, c("X", "Y"))

    testthat::expect_warning(
        out_ff <- cytobench:::ff_update_names(
            flowframe = ff,
            markermap = mk_test_markermap()
        ),
        "Could not extract spillover matrices: No spillover matrix stored in that sample"
    )
    expect_identical(flowCore::colnames(out_ff), c("X", "Y"))
    expect_identical(unname(flowCore::markernames(out_ff)), c("CD3", "CD4"))
    expect_identical(names(flowCore::markernames(out_ff)), c("X", "Y"))
})

test_that("get_new_parameternames maps named spillover columns directly", {
    ff <- mk_test_ff()
    spillmat <- matrix(c(1, 0.1, 0.1, 1), ncol = 2)
    colnames(spillmat) <- c("A", "B")

    out <- cytobench:::get_new_parameternames(
        ff = ff,
        spillmat = spillmat,
        new_channel_names = c("X", "Y"),
        markermap_names = c("A", "B")
    )

    expect_identical(out, c("X", "Y"))
    testthat::expect_warning(
        out_ff <- cytobench:::ff_update_names(
            flowframe = ff,
            markermap = mk_test_markermap()
        ),
        "Could not extract spillover matrices: No spillover matrix stored in that sample"
    )
    expect_identical(flowCore::colnames(out_ff), c("X", "Y"))
    expect_identical(unname(flowCore::markernames(out_ff)), c("CD3", "CD4"))
    expect_identical(names(flowCore::markernames(out_ff)), c("X", "Y"))
})


test_that("ff_update_names with named spillovermat", {
    ff <- mk_test_ff_spill_named()
    so_old <- flowCore::spillover(ff)
    testthat::expect_no_warning(
        out_ff <- cytobench:::ff_update_names(flowframe = ff, markermap = mk_test_markermap())
    )
    so_new <- flowCore::spillover(out_ff)
    expect_identical(unname(so_old[[3]]), unname(so_new[[3]]))
    expect_identical(flowCore::colnames(out_ff), c("X", "Y"))
    expect_identical(unname(flowCore::markernames(out_ff)), c("CD3", "CD4"))
    expect_identical(names(flowCore::markernames(out_ff)), c("X", "Y"))

    expect_equal(colnames(so_old[[3]]), c("A", "B"))
    expect_equal(colnames(so_new[[3]]), c("X", "Y"))
})


test_that("ff_update_names with number-named spillovermat", {
    ff <- mk_test_ff_spill_named()
    colnames(flowCore::keyword(ff)[["$SPILLOVER"]]) <- c("1", "2")
    so_old <- flowCore::spillover(ff)
    testthat::expect_no_warning(
        out_ff <- cytobench:::ff_update_names(flowframe = ff, markermap = mk_test_markermap())
    )
    so_new <- flowCore::spillover(out_ff)
    expect_identical(unname(so_old[[3]]), unname(so_new[[3]]))
    expect_identical(flowCore::colnames(out_ff), c("X", "Y"))
    expect_identical(unname(flowCore::markernames(out_ff)), c("CD3", "CD4"))
    expect_identical(names(flowCore::markernames(out_ff)), c("X", "Y"))

    expect_equal(colnames(so_old[[3]]), c("1", "2"))
    expect_equal(colnames(so_new[[3]]), c("X", "Y"))
})
