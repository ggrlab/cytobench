#' Simulate a FlowFrame with Normally Distributed Data
#'
#' Creates a synthetic `flowFrame` with specified column names and number of cells.
#' All values are drawn from a standard normal distribution.
#' Intended for testing or demonstration purposes only.
#'
#' @param ncells Integer. Number of simulated cells (rows).
#' @param columns Character vector of column (channel) names. Defaults to typical flow cytometry channels.
#' @param flowcore Logical. If `TRUE`, returns a `flowFrame`; if `FALSE`, returns a `data.table`.
#'
#' @return A `flowFrame` object containing random normally distributed data.
#'
#' @examples
#' ff <- simulate_ff(100)
#' flowCore::exprs(ff)[1:5, ]
simulate_ff <- function(
    ncells,
    columns = c(
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
    ),
    flowcore = TRUE) {
    data <- matrix(rnorm(ncells * length(columns)), ncol = length(columns))
    colnames(data) <- columns
    if (flowcore) {
        flowCore::flowFrame(data)
    } else {
        data.table::as.data.table(data)
    }
}


#' Simulate a FlowSet with Multiple Synthetic Samples
#'
#' Generates a `flowSet` containing `n_samples` synthetic `flowFrame`s,
#' each populated with normally distributed random values using `simulate_ff()`.
#' Intended for testing and internal development.
#'
#' @param n_samples Integer. Number of synthetic samples to generate.
#' @inheritParams simulate_ff
#'
#' @return A `flowSet` object with `n_samples` simulated `flowFrame`s.
#'
#' @examples
#' fs <- simulate_fs(3)
#' fs[[1]]
#'
#' @keywords internal
simulate_fs <- function(n_samples, flowcore = TRUE, ...) {
    tmp <- sapply(
        paste0("simsample_", seq_len(n_samples)),
        simplify = FALSE,
        function(x) simulate_ff(100, flowcore = flowcore, ...)
    )

    if (flowcore) {
        return(flowCore::flowSet(tmp))
    } else {
        return(tmp)
    }
}
