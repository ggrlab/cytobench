#' Extract Single Stain Median Fluorescence Intensity (MFI)
#'
#' This function extracts the median fluorescence intensity (MFI) from single-stain FCS files.
#'
#' @param fcs_dir A character string specifying the directory containing the FCS files. Default is "data-raw/s001".
#' @param regex_singlestain A regular expression pattern to identify single-stain FCS files. Default is "-(CD3-.*)|(none)\\.fcs".
#' @param transform
#' A function to transform the fluorescence values. Default is `function(x) { asinh(x / 1e3) }`.
#' The reported MFIs are calculated as the median of the UNtransformed values. Transformation is
#' only used to cluster the negative and positive populations.
#' @param gating_set_file A gating set file or object to apply to the FCS files. Default is NULL.
#' Usually you want to gate lymphocytes before extracting the MFIs. If a character string is provided,
#' the gating set is loaded from the file. If a (GatingSet) object is provided, it is used directly.
#'
#' @param gate_extract A specific gate to extract from the gating set. Default is NULL.
#' @return A data frame with the extracted MFIs. E.g.:
#' \preformatted{
#' # A tibble: 10 x 4
#'    feature negative positive unstained
#'    <chr>      <dbl>    <dbl>     <dbl>
#'  1 FITC-A      530.   58567.     576.
#'  2 PE-A        526.  149506.     511.
#'  3 ECD-A       453.  107584.     454.
#'  4 PC5.5-A     233.  242867.     166.
#'  5 PC7-A       176.  195307.     138.
#'  6 APC-A       136.  105036.     122.
#'  7 AF700-A     131.   28408.      97.2
#'  8 AA750-A     324.   75247.     296.
#'  9 PB-A        469.   13426.     442.
#' 10 KrO-A       442.    4473.     444.
#' }
#' @examples
#' \dontrun{
#' extracted_mfis <- extract_singlestain_mfi(
#'     "data-raw/s001",
#'     gating_set_file = "data-raw/CT_p001_s001_raw_ungated_none_Inf_navios_01-CD3-FITC-single_CD3_compensation.flowWorkspace_gatingset"
#' )
#' }
#' @export
extract_singlestain_mfi <- function(fcs_dir = "data-raw/s001",
                                    regex_singlestain = "-(CD3-.*)|(none)\\.fcs",
                                    transform = function(x) {
                                        asinh(x / 1e3)
                                    },
                                    gating_set_file = NULL,
                                    gate_extract = NULL) {
    # List all files in the specified directory that match the given regex pattern
    dir_files <- list.files(fcs_dir, full.names = TRUE, pattern = regex_singlestain)

    # Load the FCS files into a list of cytosets
    loaded_fcs <- sapply(dir_files, flowWorkspace::load_cytoset_from_fcs, simplify = FALSE)

    # Check if a gating set file is provided
    if (!all(is.null(gating_set_file))) {
        # Load the gating set from file if it's a character string, otherwise use the provided object
        if ("character" %in% class(gating_set_file)) {
            gatingset <- flowWorkspace::load_gs(gating_set_file)
        } else {
            gatingset <- gating_set_file
        }

        # Apply the gating set to each loaded FCS file
        loaded_fcs <- lapply(loaded_fcs, function(ff_x) {
            suppressMessages(
                ff_x_gated <- flowWorkspace::gh_apply_to_cs(gatingset, ff_x, compensation_source = "none")
            )

            # If no specific gate is provided, use the first actual gate (after root)
            if (is.null(gate_extract)) {
                # nolint start
                # [1] "root"
                # [2] "/FITC Control: Control Population"
                # [3] "/FITC Control: Control Population/FITC Control: FITC"
                # [4] "/FITC Control: Control Population/FITC Control: FITC/FITC Control: PE"
                # [5] "/FITC Control: Control Population/FITC Control: FITC/FITC Control: ECD"
                # [6] "/FITC Control: Control Population/FITC Control: FITC/FITC Control: PC5.5"
                # nolint end
                gate_extract <- flowWorkspace::gh_get_pop_paths(ff_x_gated)[2]
                warning("No gate_extract provided, using ", gate_extract)
            }

            # Extract the gated data
            gated_ff <- flowWorkspace::gs_pop_get_data(ff_x_gated, gate_extract) |>
                flowWorkspace::realize_view()
            return(gated_ff)
        })
    }

    # Extract median fluorescence intensities (MFIs) from the loaded FCS files
    relevant_mfis <- sapply(names(loaded_fcs), function(f_x) {
        ff_x <- flowWorkspace::cytoframe_to_flowFrame(loaded_fcs[[f_x]][[1]])
        nonempty_channel <- which(flowCore::markernames(ff_x) != "empty")

        # Check if there is more than one non-empty channel
        if (length(nonempty_channel) > 1) {
            stop("More than one non-empty channel in ", f_x)
        }

        # If no non-empty channel, it is the unstained sample
        if (length(nonempty_channel) == 0) {
            mfis <- apply(flowCore::exprs(ff_x), 2, median)
            mfis_tib <- tibble::tibble(
                "feature" = names(mfis),
                "unstained" = mfis,
            )
        } else {
            # Cluster the relevant channel into two populations and return the median of both
            values_nonempty <- flowCore::exprs(ff_x)[, names(nonempty_channel)]
            clustering <- stats::kmeans(transform(values_nonempty), centers = 2)
            mfis <- sort(tapply(values_nonempty, clustering$cluster, median))
            mfis_tib <- tibble::tibble(
                "feature" = names(nonempty_channel),
                "negative" = mfis[1],
                "positive" = mfis[2]
            )
        }

        return(mfis_tib)
    }, simplify = FALSE)

    # Combine single stainings into a single data frame
    single_stainings <- do.call(rbind, relevant_mfis[!sapply(relevant_mfis, function(x) nrow(x) > 1)])

    # Extract unstained samples
    relevant_unstained <- relevant_mfis[sapply(relevant_mfis, function(x) nrow(x) > 1)]

    # Check if there is only one unstained sample
    if (length(relevant_unstained) == 1) {
        names(relevant_unstained) <- "unstained"
    } else {
        warning("More than one unstained sample found, returning each MFI by filename")
    }

    # Join all unstained samples
    all_unstained_joint <- Reduce(dplyr::left_join, relevant_unstained)

    # Return the joined data frame of single stainings and unstained samples
    joint_df <- dplyr::left_join(
        single_stainings,
        all_unstained_joint,
        by = "feature"
    )
    # Unname each column
    return(dplyr::mutate(joint_df, across(everything(), ~unname(.))) )
}
