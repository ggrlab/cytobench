plot_markers_pairwise_base <- function(df,
                                           cofactor_namedvec,
                                           special_cofactor_list = list(),
                                           col = grDevices::rgb(0, 0, 0, .2),
                                           transform_fun = asinh,
                                           verbose = FALSE,
                                           ...) {
    # All pairwise combinations of markers
    all_combos <- utils::combn(colnames(df), 2, simplify = FALSE)
    all_combos_str <- lapply(all_combos, paste0, collapse = "_")

    df_mat <- as.matrix(df)
    if (!missing(cofactor_namedvec)) {
        df_mat[, names(cofactor_namedvec)] <- df_mat[, names(cofactor_namedvec)] %*% (diag(1 / cofactor_namedvec))
        df_mat <- transform_fun(df_mat)
    }

    markernames <- colnames(df_mat)
    par(mfrow = c(ncol(df_mat), ncol(df_mat)) - 1, mar = c(2, 2, 0, 0), oma = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
    # Loop through each pair and produce a plot
    for (marker_y in markernames[-1]) {
        for (marker_x in markernames) {
            if (verbose) {
                cat(marker_x, "vs", marker_y, ": ")
            }
            xy_str <- paste0(marker_x, "_", marker_y)
            if (marker_x == marker_y) {
                if (verbose) {
                    cat("diagonal\n")
                }
                # next
            } else if (!xy_str %in% all_combos_str) {
                if (verbose) {
                    cat("empty\n")
                }
                plot.new()
            } else {
                if (verbose) {
                    cat("plotting\n")
                }
                if (xy_str %in% names(special_cofactor_list)) {
                    vals_x <- transform_fun(df[[marker_x]] / special_cofactor_list[[xy_str]][1])
                    vals_y <- transform_fun(df[[marker_y]] / special_cofactor_list[[xy_str]][2])
                } else {
                    vals_x <- df_mat[, marker_x]
                    vals_y <- df_mat[, marker_y]
                }
                scattermore::scattermoreplot(
                    vals_x, vals_y,
                    xlab = marker_x,
                    ylab = marker_y,
                    mgp = c(0, 0, 0),
                    frame.plot = TRUE,
                    col = col,
                    ...
                )
            }
        }
    }
}
