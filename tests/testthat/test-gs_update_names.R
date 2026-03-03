devtools::load_all()

mk_test_gs <- function(channel_names = c("A", "B")) {
    exprs <- matrix(c(1, 2, 3, 4), ncol = 2)
    colnames(exprs) <- channel_names
    ff <- flowCore::flowFrame(exprs)
    tf <- tempfile(fileext = ".fcs")
    flowCore::write.FCS(ff, tf)
    cs <- flowWorkspace::load_cytoset_from_fcs(tf)
    flowWorkspace::GatingSet(cs)
}

mk_test_markermap_gs <- function(name = c("A", "B"),
                                 name_dataset1 = NULL,
                                 colnames_new = c("X", "Y")) {
    out <- data.frame(
        name = name,
        colnames = colnames_new,
        description_channel.flurochrome.marker.awh = c("CD3", "CD4"),
        stringsAsFactors = FALSE
    )
    if (!is.null(name_dataset1)) out$name_dataset1 <- name_dataset1
    out
}

mk_test_ff <- function(channel_names = c("A", "B")) {
    exprs <- matrix(c(1, 2, 3, 4), ncol = 2)
    colnames(exprs) <- channel_names
    flowCore::flowFrame(exprs)
}

test_that("gs_update_names renames channels and marker names via name", {
    gs <- mk_test_gs(c("A", "B"))
    map <- mk_test_markermap_gs()

    out <- gs_update_names(gs = gs, markermap = map)

    expect_s4_class(out, "GatingSet")
    expect_identical(flowCore::colnames(out), c("X", "Y"))
    expect_identical(unname(flowCore::markernames(out)), c("CD3", "CD4"))
    expect_identical(names(flowCore::markernames(out)), c("X", "Y"))
})

test_that("gs_update_names supports custom mapping columns", {
    gs <- mk_test_gs(c("A1", "B1"))
    map <- data.frame(
        old_custom = c("A1", "B1"),
        new_custom = c("X", "Y"),
        desc_custom = c("CD3", "CD4"),
        stringsAsFactors = FALSE
    )

    out <- gs_update_names(
        gs = gs,
        markermap = map,
        map_oldname = "old_custom",
        map_newname = "new_custom",
        map_description = "desc_custom"
    )

    expect_identical(flowCore::colnames(out), c("X", "Y"))
    expect_identical(unname(flowCore::markernames(out)), c("CD3", "CD4"))
})

test_that("gs_update_names errors on missing map columns", {
    gs <- mk_test_gs(c("A", "B"))
    bad_map <- data.frame(name = c("A", "B"), colnames = c("X", "Y"))

    expect_error(
        gs_update_names(gs = gs, markermap = bad_map),
        "missing required columns"
    )
})

test_that("gs_update_names errors when channels are not mappable", {
    gs <- mk_test_gs(c("A", "B"))
    map <- mk_test_markermap_gs(name = c("C", "D"), name_dataset1 = c("E", "F"))

    expect_message(
        expect_error(
            gs_update_names(gs = gs, markermap = map),
            "missing from marker map names"
        ),
        "Mapping of old_channel_names failed using `name`"
    )
})

test_that("gs_update_names reports useful channel mismatch details", {
    gs <- mk_test_gs(c("A", "B"))
    map <- mk_test_markermap_gs(name = c("A", "C"))

    expect_message(
        expect_error(
            gs_update_names(gs = gs, markermap = map),
            "missing from marker map names"
        ),
        "Missing: \\[B\\]"
    )
    expect_message(
        expect_error(
            gs_update_names(gs = gs, markermap = map),
            "missing from marker map names"
        ),
        "Available: \\[A, C\\]"
    )
})

test_that("gs_update_names errors when configured old-name column is absent", {
    gs <- mk_test_gs(c("A", "B"))
    map <- mk_test_markermap_gs()

    expect_error(
        gs_update_names(gs = gs, markermap = map, map_oldname = "does_not_exist"),
        "missing required columns"
    )
})

test_that("gs_update_names errors on duplicated target channel names", {
    gs <- mk_test_gs(c("A", "B"))
    map <- mk_test_markermap_gs(colnames_new = c("X", "X"))

    expect_error(
        gs_update_names(gs = gs, markermap = map),
        "contain NA or duplicates"
    )
})

test_that("check_updatednames_gs_ff validates application on a flowFrame", {
    ff <- mk_test_ff(c("A", "B"))
    gs <- mk_test_gs(c("A", "B"))
    gs <- gs_update_names(gs = gs, markermap = mk_test_markermap_gs())
    ff <- ff_update_names(flowframe = ff, markermap = mk_test_markermap_gs())

    expect_invisible(is_updatednames_gs_ff(gs = gs, ff = ff))
    expect_true(is_updatednames_gs_ff(gs = gs, ff = ff))
})
