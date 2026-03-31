test_that("Write an FCS file", {
    ff <- simulate_ff()
    tmp <- write_memory_FCS(ff)
    reread_ff <- flowCore::read.FCS(tmp)
    testthat::expect_false(identical(reread_ff, ff)) # they SHOULD be different objects!
    testthat::expect_equal(nrow(reread_ff), nrow(ff))
    testthat::expect_true(max(abs(flowCore::exprs(ff) - flowCore::exprs(reread_ff))) < 1e-8)
})

test_that("write_memory_FCS uses /dev/shm when available", {
    skip_if_not(dir.exists("/dev/shm"), message = "/dev/shm not available")
    ff <- simulate_ff()
    tmp <- write_memory_FCS(ff)
    testthat::expect_true(startsWith(tmp, "/dev/shm/"))
})

test_that("write_memory_FCS cleans up when envir exits", {
    ff <- simulate_ff()
    path <- local({
        tmp <- write_memory_FCS(ff, envir = environment())
        testthat::expect_true(file.exists(tmp))
        tmp
    })
    testthat::expect_false(file.exists(path))
})

test_that("write_memory_FCS respects local envir lifetime", {
    ff <- simulate_ff()
    path <- local({
        tmp <- write_memory_FCS(ff, envir = environment())
        testthat::expect_true(file.exists(tmp))
        tmp
    })
    # File does NOT survive the local{} block because it's tied to global env
    testthat::expect_false(file.exists(path))
})

test_that("write_memory_FCS keeps file alive with global envir", {
    ff <- simulate_ff()
    path <- local({
        tmp <- write_memory_FCS(ff, envir = rlang::global_env())
        testthat::expect_true(file.exists(tmp))
        tmp
    })
    # File survives the local{} block because it's tied to global env
    testthat::expect_true(file.exists(path))
    # Manual cleanup
    unlink(path, force = TRUE)
})

test_that("write_memory_FCS keeps file alive with global envir: THIS IS DEFAULT!", {
    ff <- simulate_ff()
    path <- local({
        tmp <- write_memory_FCS(ff)
        testthat::expect_true(file.exists(tmp))
        tmp
    })
    # File survives the local{} block because it's tied to global env
    testthat::expect_true(file.exists(path))
    # Manual cleanup
    unlink(path, force = TRUE)
})

test_that("write_memory_FCS lifetime follows envir, not call site", {
    ff <- simulate_ff()

    wrapper <- function() {
        helper <- function() {
            write_memory_FCS(ff, envir = parent.frame())
        }
        path <- helper()
        # File survives helper() because it's tied to wrapper()
        testthat::expect_true(file.exists(path))
        path
    }

    path <- wrapper()
    # File is gone because wrapper()'s frame exited
    testthat::expect_false(file.exists(path))
})
