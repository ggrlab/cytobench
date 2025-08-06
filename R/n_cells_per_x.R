#' Calculate the number of cells per specified grouping
#'
#' Computes the number of cells for each combination of ID columns (e.g., sample)
#' and grouping variables (e.g., cluster, metaCluster), returning wide-format tables
#' for each grouping. Ensures that all levels are represented even if missing in a sample.
#'
#' @param dt A `data.table` or `data.frame` containing cell-level data with grouping columns.
#' @param which Character vector. The grouping variables to summarize (e.g., \code{c("cluster", "metaCluster")}).
#' @param id_cols Character vector. Identifier columns to group by (e.g., \code{"sample"}).
#' @param which_levels Optional named list of levels per grouping variable.
#' If \code{NULL}, the levels are inferred from the data by:
#' 1. Converting the unique values of the column to numeric.
#' 2. Taking the maximum value (max_level_numeric)
#' 3. Setting the levels to 1:max_level_numeric
#'
#' @return A named list. Each entry corresponds to one grouping variable from \code{which},
#' containing a wide-format data frame with counts per level (e.g., one row per sample, one column per cluster).
#'

#' @return A data frame in wide format with the counts of cells per specified grouping.
#'
#' @examples
#' # Example usage:
#' dt <- data.table::data.table(
#'     sample = as.character(rep(1:3, each = 4)),
#'     cluster = rep(1:2, times = 6),
#'     metaCluster = rep(1:4, times = 3)
#' )
#' n_cells_per_x(dt)
#' @export
#' @keywords flowsom
n_cells_per_x <- function(dt,
                          which = c("cluster", "metaCluster"),
                          id_cols = c("sample"),
                          which_levels = NULL) {
    # Ensure input is a data.table
    if (!data.table::is.data.table(dt)) {
        dt <- data.table::as.data.table(dt)
    }

    # For each grouping variable, compute counts and return wide-format
    sapply(which, function(x) {
        grouping_columns <- c(id_cols, x)

        # Count number of cells per ID/group
        tmp <- dt[, .N, by = grouping_columns]

        # If no levels given, infer them from the data
        if (is.null(which_levels[[x]])) {
            max_level_numeric <- max(as.numeric(unique(dt[[x]])))
            which_levels[[x]] <- 1:max_level_numeric
        } else if (length(which_levels[[x]]) == 1) {
            which_levels[[x]] <- 1:which_levels[[x]]
        }

        # Convert group column to ordered factor with level labels
        tmp[[x]] <- factor(tmp[[x]],
            levels = sort(as.numeric(which_levels[[x]]))
        )
        levels(tmp[[x]]) <- paste0(x, "_", levels(tmp[[x]]))

        # Add dummy row to enforce presence of all factor levels
        dummy_row <- data.frame("_____REMOVEME_____", levels(tmp[[x]]), 0)
        colnames(dummy_row) <- colnames(tmp)

        # Ensure full level coverage via dummy row; drop it afterwards
        tmp_wide <- dplyr::add_row(
            tmp |> tibble::as_tibble(),
            # Add a REMOVEME sample just to make sure that all clusters are present in the output
            dummy_row
        ) |>
            tidyr::pivot_wider(
                names_from = tidyr::all_of(x),
                values_from = "N",
                values_fill = 0
            ) |>
            dplyr::filter(sample != "_____REMOVEME_____")

        # Return only ID columns and count columns (reordered)
        tmp_wide[, c(id_cols, levels(tmp[[x]]))]
    }, USE.NAMES = TRUE, simplify = FALSE)
}
