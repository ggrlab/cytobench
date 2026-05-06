#' Extract per-cell gate membership (and optionally expression values) from a GatingHierarchy
#'
#' For a `GatingHierarchy`, retrieves the logical
#' gate-membership vector for one or more populations and binds them into a
#' single `data.table`. Optionally appends the per-cell expression matrix so
#' that each row corresponds to one event annotated with both its measured
#' channel values and its membership in each requested gate.
#'
#' Gate membership is obtained via [flowWorkspace::gh_pop_get_indices()],
#' which returns a logical vector of length equal to the number of events
#' in the sample.
#'
#' @param gatinghierarchy A [`flowWorkspace::GatingHierarchy`] object.
#' @param gatenames Character vector of gate paths (e.g. `"/root/Lymph/CD4"`),
#'   integer indices into the population paths returned by
#'   [flowWorkspace::gs_get_pop_paths()], or a mix of the two. If `NULL`
#'   (the default), all population paths in the gating hierarchy are used.
#' @param return_exprs Logical. If `TRUE`, the per-event expression matrix
#'   is column-bound to the gate-membership columns. If `FALSE` (default),
#'   only the gate-membership columns are returned.
#'
#' @return A [`data.table::data.table`] with one row per event. Columns are
#'   the resolved gate paths (logical), optionally prefixed by the channels
#'   from the expression matrix when `return_exprs = TRUE`.
#'
#' @examples
#' \dontrun{
#' # All gates, membership only
#' membership <- assign_gatingset_cells(gs)
#'
#' # Specific gates by path and index, with expression values
#' full <- assign_gatingset_cells(
#'     gs,
#'     gatenames = list("/root/Lymph/CD4", 5L),
#'     return_exprs = TRUE
#' )
#' }
#'
#' @export
#' @keywords cytometry
assign_gatingset_cells <- function(gatinghierarchy, gatenames = NULL, return_exprs = FALSE) {
    # Resolve numeric gate identifiers to gate paths.
    pop_paths <- flowWorkspace::gs_get_pop_paths(gatinghierarchy)
    if (is.null(gatenames)) {
        gatenames <- pop_paths
    }
    gatenames <- sapply(
        gatenames, function(gatename) {
            if (is.numeric(gatename)) {
                pop_paths[[gatename]]
            } else {
                gatename
            }
        }
    )
    final_extracted <- sapply(gatenames, simplify = FALSE, function(gatename) {
        data.table::as.data.table(flowWorkspace::gh_pop_get_indices(gatinghierarchy, gatename))
    }) |> data.table::cbindlist()
    colnames(final_extracted) <- gatenames

    if (!return_exprs) {
        return(final_extracted)
    }
    exprs_list <- lapply(seq_along(gatinghierarchy), function(i) {
        flowWorkspace::gs_get_cytoframe(gatinghierarchy, i) |>
            flowCore::exprs() |>
            data.table::as.data.table()
    })
    exprs_concatenated <- data.table::rbindlist(exprs_list)

    return(cbind(exprs_concatenated, final_extracted))
}
