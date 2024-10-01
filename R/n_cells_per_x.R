#' Calculate the number of cells per specified grouping
#'
#' This function calculates the number of cells per specified grouping (e.g., cluster or metaCluster)
#' and returns a wide-format data frame with the counts.
#'
#' @param dt A data.table containing the data.
#' @param which A character vector specifying the columns to group by. Default is c("cluster", "metaCluster").
#' @param id_cols A character vector specifying the ID columns. Default is c("sample").
#' @param which_levels A list specifying the levels for each grouping column. Default is NULL.
#' If NULL, the levels are inferred from the data by:
#' 1. Converting the unique values of the column to numeric.
#' 2. Taking the maximum value (max_level_numeric)
#' 3. Setting the levels to 1:max_level_numeric
#'
#' @return A data frame in wide format with the counts of cells per specified grouping.
#' @importFrom dplyr add_row filter
#' @importFrom tidyr pivot_wider all_of
#' @export
#'
#' @examples
#' # Example usage:
#' dt <- data.table::data.table(
#'     sample = as.character(rep(1:3, each = 4)),
#'     cluster = rep(1:2, times = 6),
#'     metaCluster = rep(1:4, times = 3)
#' )
#' n_cells_per_x(dt)
n_cells_per_x <- function(dt,
                          which = c("cluster", "metaCluster"),
                          id_cols = c("sample"),
                          which_levels = NULL) {
    if (!data.table::is.data.table(dt)) {
        data.table::as.data.table(dt)
    }
    sapply(which, function(x) {
        grouping_columns <- c(id_cols, x)
        tmp <- dt[, .N, by = grouping_columns]
        if (is.null(which_levels[[x]])) {
            max_level_numeric <- max(as.numeric(unique(dt[[x]])))
            which_levels[[x]] <- 1:max_level_numeric
        } else if (length(which_levels[[x]]) == 1) {
            which_levels[[x]] <- 1:which_levels[[x]]
        }
        tmp[[x]] <- factor(tmp[[x]], levels = sort(as.numeric(which_levels[[x]])))
        levels(tmp[[x]]) <- paste0(x, "_", levels(tmp[[x]]))

        removin_sample <- data.frame(
            "_____REMOVEME_____", levels(tmp[[x]]), 0
        )
        # ensure that the column names are the same
        colnames(removin_sample) <- colnames(tmp)
        tmp_wide <- dplyr::add_row(
            tmp |> tibble::as_tibble(),
            # Add a REMOVEME sample just to make sure that all clusters are present in the output
            removin_sample
        ) |>
            tidyr::pivot_wider(
                names_from = tidyr::all_of(x),
                values_from = "N",
                values_fill = 0
            ) |>
            # After the pivot_wider, the REMOVEME sample is not needed anymore
            dplyr::filter(sample != "_____REMOVEME_____")
        # The following mainly resorts the cluster_I columns
        tmp_wide[, c(id_cols, levels(tmp[[x]]))]
    }, USE.NAMES = TRUE, simplify = FALSE)
}
