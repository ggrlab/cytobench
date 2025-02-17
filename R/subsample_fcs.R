#' Subsample FCS File
#'
#' This function reads an FCS file, subsamples it to a specified number of cells,
#' and writes the subsampled data to a new file.
#'
#' @param file_x Character. Path to the input FCS file.
#' @param n_cells Integer. Number of cells to subsample. Default is 10,000.
#' @param seed Integer. Seed for random number generation to ensure reproducibility. Default is 427764.
#' @param outdir Character. Directory where the subsampled FCS file will be saved. Default is "res/a01/02_subsampled/random_<n_cells>".
#' @param indir Character. Directory where the original FCS file is located. Default is "res/a04/01_gated".
#' @param verbose Logical. If TRUE, prints progress messages. Default is TRUE.
#'
#' @return None. The function writes the subsampled FCS file to the specified output directory.
#'
#' @examples
#' \dontrun{
#' subsample_fcs("path/to/input.fcs", n_cells = 5000, seed = 12345, outdir = "output/directory", indir = "input/directory", verbose = TRUE)
#' }
#'
#' @export
subsample_fcs <- function(file_x,
                          n_cells = 10000,
                          seed = 427764,
                          outdir = paste0("res/a01/02_subsampled/random_", n_cells),
                          indir = "res/a04/01_gated",
                          verbose = TRUE) {
    if (verbose) {
        cat("Processing", file_x, "...\n")
    }
    ff <- flowCore::read.FCS(file_x)
    ff_subsampled <- subsample_ff(ff, seed = seed, n_cells = n_cells)

    file_new <- sub(indir, "", file_x)
    file_new <- sub("_none_Inf", paste0("_random_", n_cells), file_new)
    path_new <- file.path(outdir, file_new)
    dir.create(dirname(path_new), recursive = TRUE, showWarnings = FALSE)
    flowCore::write.FCS(ff_subsampled, path_new)
    if (verbose) {
        cat("   Wrote to", path_new, "\n")
    }
}


#' Subsample or Upsample a FlowFrame
#'
#' This function subsamples or upsamples a given FlowFrame object to a specified number of cells.
#' If the number of cells in the FlowFrame is less than the specified number, the function
#' replicate the cells until the specified number is reached.
#' E.g. given 100 cells and n_cells = 5000, the function will replicate each cell of the
#' 100 cells exactly 50 times.
#'
#' @param flowframe A FlowFrame object containing the data to be subsampled or upsampled.
#' @param n_cells An integer specifying the number of cells to sample. Default is 10,000.
#' @param seed An integer specifying the seed for random number generation. Default is 427764.
#'
#' @return A FlowFrame object containing the subsampled or upsampled data.
#'
#' @examples
#' # Assuming `ff` is a FlowFrame object with flow cytometry data
#' subsampled_ff <- subsample_ff(ff, n_cells = 5000, seed = 12345)
#'
#' @export
subsample_ff <- function(flowframe,
                         n_cells = 10000,
                         seed = 427764) {
    # Note that if you select more cells than are available, the function will upsample.
    # But the first indices will be the same as the original flowframe. (up to the number
    # of cells in the original flowframe)
    seq_ff <- seq_len(nrow(flowframe))
    set.seed(seed)
    # Sample n_cells from the flowframe.
    if (nrow(flowframe) >= n_cells) {
        subsampled_indices <- sample(seq_ff, n_cells, replace = FALSE)
    } else {
        # Then we have to UPsample.
        # Cell selection similar to torch_geometric's fixed_points
        # https://pytorch-geometric.readthedocs.io/en/latest/_modules/torch_geometric/transforms/fixed_points.html

        # 1. Take all cells
        index_permutations <- seq_ff
        # 2. How often do we have to replicate the cells to reach n_cells?
        n_replications <- ceiling(n_cells / nrow(flowframe))

        index_permutations <- cbind(
            index_permutations,
            replicate(
                n_replications,
                # For each replicate, we permute the indices WITHOUT replacement
                sample(seq_ff, replace = FALSE)
            )
        )
        # 3. Flatten the matrix
        subsampled_indices <- unlist(index_permutations)[1:n_cells]
        if (!all(subsampled_indices[seq_ff] == seq_ff)) {
            stop("Subsampling failed - The initial indices are not in the first nrow(ff) indices.")
        }
    }

    # Sub/Up-sample the flowframe
    ff_subsampled <- flowframe[subsampled_indices, ]
    return(ff_subsampled)
}
