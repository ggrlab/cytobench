#' Gate Cells Based on Compensated FCS and Gating Set Files
#'
#' This function processes compensated FCS files and their corresponding gating set files to gate cells based on specified criteria.
#'
#' @param dir_fcs_compensated Character. Directory containing compensated FCS files. Note that any compensation in these files is NOT applied.
#' @param dir_gatingsets_compensated Character. Directory containing compensated gating set files.
#' The gating set files are matched to each FCS file based on the sID (sampleID) from FCS and Gating *directory* name.
#' @param pattern_files Character. Pattern to match FCS files. Default is "panel".
#' @param outdir Character. Output directory to save results. Default is NA, then no output is saved.
#' @param verbose Logical. If TRUE, prints progress messages. Default is TRUE.
#' @param gatename Character or Numeric. Name or index of the gate to extract data from. Default is "/Singlets/CD45+/CD3+".
#' @param filename_ungated_into Character. String to replace "_ungated_" in the filename. Default is "_cd3_".
#' @param exclude_pops Character vector. Populations to exclude from the final counts. Default is "/Singlets/Lymphocytes", "/Singlets/Monocytes", "/Singlets/Granulocytes".
#'
#' @return A data.table containing the counts and percentages of gated populations. Usually you would use "outdir" as a side effect to save the data.
#' @export
#'
#' @examples
#' \dontrun{
#' gate_cells(
#'     dir_fcs_compensated = "path/to/compensated/fcs",
#'     dir_gatingsets_compensated = "path/to/gatingsets",
#'     outdir = "path/to/output"
#' )
#' }
gate_cells <- function(flowset,
                       gatingset,
                       gatename = "/Singlets/CD45+/CD3+",
                       verbose = TRUE,
                       inplace = FALSE) {
    global_gating <- NULL
    if (length(gatingset) <= 1) {
        global_gating <- gatingset
    }

    if (is.null(global_gating)) {
        # 3. Check if the panel files and gating files match
        # 3.1 Get the sample IDs from the directory names of the panel files
        fcs_panel_ids <- flowCore::sampleNames(flowset)
        gates_panel_ids <- names(gatingset)
        if (anyDuplicated(gates_panel_ids) != 0) {
            stop("There are duplicated sample IDs in the gating files, cannot assign them uniquely to the to-be-gated FCS files")
        }

        # 3.2 Check if all panel files have a corresponding gating file
        if (!all(fcs_panel_ids %in% gates_panel_ids)) {
            # In the case of CT and MX data, each sample was measured on NAVIOS and Cytoflex,
            # thus it could be that the gating file of NAVIOS is missing.
            # However, this would then stop with an error in the following for-loop
            stop("Not all panel files have a corresponding gating file")
        }
    }
    if(!inplace){
        flowset <- flowCore::flowSet(flowCore::flowSet_to_list(flowset))
    }
    # Initialize an empty variable to store complete counts
    counts_complete <- NA
    # Loop through each compensated FCS panel file
    for (fcs_i in seq_along(flowset)) {
        fcs_x <- flowset[[fcs_i]]
        fcs_x_sID <- flowCore::sampleNames(flowset)[fcs_i]
        if (verbose) {
            cat("Processing ", fcs_x_sID, "\n")
        }
        if (is.null(global_gating)) {
            gating_x <- gatingset[fcs_x_sID]
            if (length(gating_x) != 1) {
                stop("No unique gate file found for ", fcs_x_sID, " in ", paste0(names(gatingset), collapse = ", "))
            }
        } else {
            gating_x <- global_gating
        }
        # Write the current fcs into a temporary directory
        ff_path_tmp <- write_memory_FCS(fcs_x)
        # Load the cytoset from the FCS file
        cs <- flowWorkspace::load_cytoset_from_fcs(ff_path_tmp)
        # Apply the gating set to the cytoset and suppress messages
        applied_gates <- suppressMessages(flowWorkspace::gh_apply_to_cs(
            gating_x,
            cs,
            # The following is NOT necessary for my "compensated.manual" fcs files
            # because their compensation matrix is 1) applied and 2) set to 0,
            # thus applying or not applying the compensation matrix here does not change
            # the data.
            compensation_source = "none"
        ))
        suppressMessages(flowWorkspace::recompute(applied_gates))
        # Get the counts and percentages for each population
        counts_freqs <- lapply(c("count", "percent"), function(statistic) {
            flowWorkspace::gs_pop_get_stats(
                applied_gates,
                type = statistic
            )
        })
        # Merge the counts and percentages into a single data table
        counts_dt <- dplyr::left_join(
            counts_freqs[[1]],
            counts_freqs[[2]],
            by = c("sample", "pop")
        )

        # Rename the percentage column
        colnames(counts_dt)[4] <- "percent_from_parent"
        # Extract the population MFIs
        # This only extracts for the columns in flowCore::markernames()
        pop_mfis <- flowWorkspace::gh_pop_get_stats(applied_gates, type = flowWorkspace::pop.MFI)
        counts_dt_MFI <- dplyr::left_join(counts_dt, pop_mfis, by = c("pop"))
        counts_dt_MFI[["sample"]] <- fcs_x_sID

        # Combine the counts with the complete counts
        if (all(is.na(counts_complete))) {
            counts_complete <- counts_dt_MFI
        } else {
            counts_complete <- dplyr::bind_rows(counts_complete, counts_dt_MFI)
        }

        # If gatename is numeric, get the corresponding population path
        if (is.numeric(gatename)) {
            gatename <- flowWorkspace::gs_get_pop_paths(gating_x)[[gatename]]
        }
        # Extract the CD3+ population data
        extracted <- flowWorkspace::gh_pop_get_data(applied_gates, gatename)
        extracted <- flowWorkspace::realize_view(extracted)
        extracted_ff <- flowWorkspace::cytoframe_to_flowFrame(extracted)

        flowset[[fcs_i]] <- extracted_ff
    }

    return(list(
        counts = counts_complete,
        flowset_gated = flowset
    ))
}
