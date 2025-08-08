testthat::test_that("Fast sweep replacement", {
    n <- 1e2
    n_colums <- 10
    n_rows <- n / n_colums
    m <- matrix(rnorm(n_colums), nrow = n_rows, ncol = n_colums)
    stats <- rnorm(n_colums)
    testthat::expect_equal(
        sweep(m, 2, stats, "-"),
        sweep_matrix_col(m, stats, "-")
    )
    testthat::expect_equal(
        sweep(m, 2, stats, "+"),
        sweep_matrix_col(m, stats, "+")
    )
    testthat::expect_equal(
        sweep(m, 2, stats, "*"),
        sweep_matrix_col(m, stats, "*")
    )
    testthat::expect_equal(
        sweep(m, 2, stats, "^"),
        sweep_matrix_col(m, stats, "^")
    )


    n <- 1e5
    n_colums <- 1e4
    n_rows <- n / n_colums
    m <- matrix(rnorm(n_colums), nrow = n_rows, ncol = n_colums)
    stats <- rnorm(n_colums)
    # microbenchmark::microbenchmark(
    #     v1 = sweep(m, 2, stats, "-"),
    #     v2 = sweep_matrix_col(m, stats, "-"),
    #     times = 100,
    #     check = "equal"
    # )


    n <- 1e5
    n_colums <- 10
    n_rows <- n / n_colums
    m <- matrix(rnorm(n_colums), nrow = n_rows, ncol = n_colums)
    stats <- rnorm(n_colums)
    # microbenchmark::microbenchmark(
    #     v1 = sweep(m, 2, stats, "-"),
    #     v2 = sweep_matrix_col(m, stats, "-"),
    #     times = 100,
    #     check = "equal"
    # )
})
