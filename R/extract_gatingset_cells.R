#' Extract Gated Population Data From a `GatingSet`
#'
#' This function retrieves one or more gated populations from a `GatingSet`
#' and returns them either as `flowFrame` objects (`flowset` option),
#' as `cytoset` views (`cytoset` option), or returns the original
#' `GatingSet` (`gatinghierarchy` option).
#'
#' Internally, it also computes population-level counts, percentages,
#' and MFIs for the provided `gatingset`. These statistics are currently
#' not returned by the function but are kept for compatibility with related
#' gating workflows in this package.
#'
#' @param gatingset A `flowWorkspace::GatingSet` object.
#' @param gatenames Character or numeric vector identifying populations
#'   to extract. Character values are interpreted as gate paths.
#'   Numeric values are interpreted as indices into
#'   `flowWorkspace::gs_get_pop_paths(gatingset)`.
#' @param return_x Character vector specifying the output type.
#'   Supported values are `"flowset"`, `"cytoset"`, and
#'   `"gatinghierarchy"`. Only the first element is used.
#'
#' @return
#' - If `return_x[1] == "gatinghierarchy"`: returns `gatingset` unchanged.
#' - Otherwise: returns a named list where each element corresponds to one
#'   entry in `gatenames`:
#'   - `"cytoset"`: each element is a realized `cytoset` view.
#'   - `"flowset"`: each element is a `flowFrame` converted from the
#'     realized `cytoset` view.
#'
#' @details
#' MFI extraction is first attempted with `inverse.transform = TRUE`.
#' If that fails, the function retries with `inverse.transform = FALSE`
#' and emits a warning.
#'
#' @export
#' @keywords cytometry
#' @examples
#' \dontrun{
#' gs <- flowWorkspace::load_gating_set("path/to/gating_set")
#'
#' # Extract root and a specific population as flowFrames
#' gated_extract_cells(
#'     gatingset = gs,
#'     gatenames = c("root", "/Singlets/CD45+/CD3+"),
#'     return_x = "flowset"
#' )
#' }
extract_gatingset_cells <- function(gatingset, gatenames = c("root"), return_x = c("flowset", "cytoset", "gatinghierarchy")) {
    # Retrieve per-population counts and percentages.
    counts_freqs <- lapply(c("count", "percent"), function(statistic) {
        flowWorkspace::gs_pop_get_stats(
            gatingset,
            type = statistic
        )
    })

    # Merge count and percent tables by sample and population.
    counts_dt <- dplyr::left_join(
        counts_freqs[[1]],
        counts_freqs[[2]],
        by = c("sample", "pop")
    )

    # Standardize percentage column naming.
    colnames(counts_dt)[4] <- "percent_from_parent"

    # Extract population MFIs for marker columns.
    pop_mfis <- tryCatch(
        {
            flowWorkspace::gh_pop_get_stats(
                gatingset,
                type = flowWorkspace::pop.MFI,
                inverse.transform = TRUE
            )
        },
        error = function(e) {
            warning(
                "Error due to inverse.transform in MFI extraction  : ", e$message
            )
            flowWorkspace::gh_pop_get_stats(
                gatingset,
                type = flowWorkspace::pop.MFI,
                inverse.transform = FALSE
            )
        }
    )
    counts_dt_MFI <- dplyr::left_join(counts_dt, pop_mfis, by = "pop")


    # Resolve numeric gate identifiers to gate paths.
    pop_paths <- flowWorkspace::gs_get_pop_paths(gatingset)
    gatenames <- sapply(gatenames, function(gatename) {
        if (is.numeric(gatename)) {
            pop_paths[[gatename]]
        } else {
            gatename
        }
    })


    # Optionally return gating hierarchy unchanged.
    if (return_x[1] == "gatinghierarchy") {
        final_extracted <- gatingset
    } else {
        # Extract requested populations and coerce output format if needed.
        final_extracted <- sapply(gatenames, simplify = FALSE, function(gatename) {
            extracted <- flowWorkspace::gh_pop_get_data(gatingset, gatename)
            extracted <- flowWorkspace::realize_view(extracted)
            if (return_x[1] == "cytoset") {
                extracted
            } else {
                flowWorkspace::cytoframe_to_flowFrame(extracted)
            }
        })
    }

    return(
        list(
            counts = counts_dt_MFI,
            extracted = final_extracted
        )
    )
}
