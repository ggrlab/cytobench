test_that("min-max scaling with two values works", {
    dt <- data.table(x = c(10, 55, 100))
    scale_column_minmax(dt, scaling_values = c(10, 100), colX = "x")
    expect_equal(dt$x, c(0, 0.5, 1))
})

test_that("scaling with one value performs centering", {
    dt <- data.table(x = c(90, 100, 110))
    scale_column_minmax(dt, scaling_values = 100, colX = "x")
    expect_equal(dt$x, c(-10, 0, 10))
})

test_that("in-place modification occurs", {
    dt <- data.table(x = 1:5)
    result <- scale_column_minmax(dt, scaling_values = c(1, 5), colX = "x")
    expect_identical(result, dt)
})

test_that("subtract_bg = FALSE disables background subtraction", {
    dt <- data.table(x = c(10, 55, 100))
    scale_column_minmax(dt, scaling_values = c(10, 100), colX = "x", subtract_bg = FALSE)
    expect_equal(dt$x, c(10, 55, 100) / 90)
})

test_that("one-value centering ignores subtract_bg", {
    dt1 <- data.table(x = c(90, 100, 110))
    dt2 <- data.table(x = c(90, 100, 110))
    scale_column_minmax(dt1, scaling_values = 100, colX = "x", subtract_bg = TRUE)
    scale_column_minmax(dt2, scaling_values = 100, colX = "x", subtract_bg = FALSE)
    expect_equal(dt1, dt2) # should be identical
})

test_that("error is thrown for more than 2 scaling values", {
    dt <- data.table(x = 1:5)
    expect_error(
        scale_column_minmax(dt, scaling_values = c(1, 2, 3), colX = "x"),
        "Expected 1 or 2 scaling values"
    )
})

test_that("scaling on column with NAs works", {
    dt <- data.table(x = c(0, NA, 1))
    scale_column_minmax(dt, scaling_values = c(0, 1), colX = "x")
    expect_equal(dt$x, c(0, NA, 1))
})

test_that("non-numeric column throws an error", {
    dt <- data.table(x = c("a", "b", "c"))
    expect_error(
        scale_column_minmax(dt, scaling_values = c(1, 2), colX = "x"),
        "non-numeric argument"
    )
})
