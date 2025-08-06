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
#' @export
#' @keywords test-helper
simulate_ff <- function(
    ncells = 250,
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
    data <- matrix(stats::rnorm(ncells * length(columns)), ncol = length(columns))
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
#' @param ... Additional arguments passed to `simulate_ff()`.
#'
#' @return A `flowSet` object with `n_samples` simulated `flowFrame`s.
#'
#' @examples
#' fs <- simulate_fs(3)
#' fs[[1]]
#'
#' @export
#' @keywords test-helper
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
#' Simulate Multi-Stain Flow Cytometry Data
#'
#' Generates a simulated `flowSet` containing single-stained, unstained, and multi-stained FCS-like samples.
#' Each single-stain sample has signal only in one marker channel; multi-stained and
#' unstained controls are also included.
#'
#' This simulated data reflects the data from our CD3 manuscript and is useful for testing.
#'
#' @param columns A character vector of marker channels to simulate (default: common 10-color panel).
#'
#' @return A `flowSet` of simulated `flowFrame` objects, with appropriate `markernames()` set
#' for downstream compatibility (e.g. in MFI extraction functions).
#'
#' @details
#' - Single-stained files are named `sample0_01-CD3-FITC.fcs`, ..., one per marker.
#' - Unstained control is named `sample0_11-none`.
#' - Panel-stained and master-mix controls are named `sample0_12-panel` and `sample0_15-MasterMix`.
#' - All simulated data is random normal noise scaled to synthetic MFI ranges.
#'
#' @examples
#' fs_sim <- simulate_cd3()
#' flowCore::sampleNames(fs_sim)
#'
#' @export
#' @keywords test-helper
simulate_cd3 <- function(columns = c(
                             "FITC-A", "PE-A", "ECD-A", "PC5.5-A", "PC7-A",
                             "APC-A", "AF700-A", "AA750-A", "PB-A", "KrO-A"
                         )) {
    fs_singlestain <- simulate_fs(length(columns) + 3, flowcore = FALSE)

    # Name samples: one for each marker, plus 3 extra controls
    samplenames_singlestain <- paste0("sample0_", sprintf("%02d-CD3-%s.fcs", seq_along(columns), sub("-A$", "", columns)))
    names(samplenames_singlestain) <- columns
    names(fs_singlestain) <- c(
        samplenames_singlestain,
        paste0("sample0_", c("11-none", "12-panel", "15-MasterMix"))
    )

    samplefun <- function(n, sample_from = c(1, 1e4)) {
        sample(sample_from, size = n, replace = TRUE) * stats::rnorm(n, mean = 10, sd = 1)
    }

    for (marker_x in columns) {
        sample_x <- samplenames_singlestain[[marker_x]]
        n_cells <- nrow(fs_singlestain[[sample_x]])
        fs_singlestain[["sample0_12-panel"]][, marker_x] <- samplefun(n_cells, sample_from = 5e3)
        fs_singlestain[[sample_x]][, marker_x] <- samplefun(n_cells)
        fs_singlestain[["sample0_15-MasterMix"]][, marker_x] <- samplefun(n_cells)
    }

    # Convert to flowFrames and set marker names
    ff_list <- lapply(fs_singlestain, function(x) {
        flowCore::flowFrame(as.matrix(x))
    })

    for (marker_x in columns) {
        sample_x <- samplenames_singlestain[[marker_x]]
        flowCore::markernames(ff_list[[sample_x]])[TRUE] <- "empty"
        flowCore::markernames(ff_list[[sample_x]])[marker_x] <- marker_x
    }
    flowCore::markernames(ff_list[["sample0_11-none"]])[TRUE] <- "empty"

    return(flowCore::flowSet(ff_list))
}
