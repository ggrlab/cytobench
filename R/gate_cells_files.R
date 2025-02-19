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
gate_cells_files <- function(dir_fcs_compensated,
                       dir_gatingsets_compensated,
                       pattern_files = "panel",
                       outdir = NA,
                       verbose = TRUE,
                       gatename = "/Singlets/CD45+/CD3+",
                       filename_ungated_into = "_cd3_",
                       exclude_pops = paste0("/Singlets/", c("Lymphocytes", "Monocytes", "Granulocytes"))) {
    # 1. List all panel files in the compensated FCS directory given pattern_files
    compensated_fcs_panel <- list.files(
        dir_fcs_compensated,
        pattern = pattern_files,
        full.names = TRUE,
        recursive = TRUE
    )

    # 2. List all gating files in the compensated gating sets directory
    # For each panel file, there should be one corresponding gating file
    compensated_gates_panel <- list.files(
        dir_gatingsets_compensated,
        pattern = ".flowWorkspace_gatingset",
        full.names = TRUE,
        recursive = TRUE,
        include.dirs = TRUE
    )

    global_gating <- NULL
    if (length(compensated_gates_panel) == 0) {
        global_gating <- tryCatch(
            flowWorkspace::load_gs(dir_gatingsets_compensated),
            error = function(e) {
                stop("No gating files found in ", dir_gatingsets_compensated)
            }
        )
    }

    if (is.null(global_gating)) {
        # 3. Check if the panel files and gating files match
        # 3.1 Get the sample IDs from the directory names of the panel files
        comp_fcs_panel_ids <- basename(dirname(compensated_fcs_panel))
        comp_gates_panel_ids <- basename(dirname(compensated_gates_panel))
        if (anyDuplicated(comp_gates_panel_ids) != 0) {
            stop("There are duplicated sample IDs in the gating files, cannot assign them uniquely to the to-be-gated FCS files")
        }
        names(compensated_gates_panel) <- comp_gates_panel_ids
        names(compensated_fcs_panel) <- comp_fcs_panel_ids

        # 3.2 Check if all panel files have a corresponding gating file
        if (!all(comp_fcs_panel_ids %in% comp_gates_panel_ids)) {
            # In the case of CT and MX data, each sample was measured on NAVIOS and Cytoflex,
            # thus it could be that the gating file of NAVIOS is missing.
            # However, this would then stop with an error in the following for-loop
            stop("Not all panel files have a corresponding gating file")
        }
    }
    # Initialize an empty variable to store complete counts
    counts_complete <- NA
    # Loop through each compensated FCS panel file
    for (fcs_i in seq_along(compensated_fcs_panel)) {
        fcs_x <- compensated_fcs_panel[fcs_i]
        if (verbose) {
            cat("Processing ", fcs_x, "\n")
        }
        fcs_x_sID <- names(compensated_fcs_panel)[fcs_i]
        if (is.null(global_gating)) {
            # The compensated_gates_panel are named without "_comp-manual" in the filename but with "_raw" instead.
            # So to find the correct gating file, we need to replace this
            gate_filename <- compensated_gates_panel[names(compensated_gates_panel) == fcs_x_sID]
            if (length(gate_filename) != 1) {
                stop("No unique gate file found for ", fcs_x, "with sID=", fcs_x_sID, " in ", paste0(names(compensated_gates_panel), collapse = ", "))
            }
            # Load the gating set from the gating file
            gating_x <- flowWorkspace::load_gs(gate_filename)
        } else {
            gating_x <- global_gating
        }
        # Load the cytoset from the FCS file
        cs <- flowWorkspace::load_cytoset_from_fcs(fcs_x)
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
        pop_mfis <- flowWorkspace::gh_pop_get_stats(applied_gates, type = flowWorkspace::pop.MFI)
        # For some reason, "sampleID" is not included in the stats
        colnames(pop_mfis)[-1] <- flowCore::colnames(applied_gates[[1]])[flowCore::colnames(applied_gates[[1]]) != "sampleID"]
        counts_dt_MFI <- dplyr::left_join(counts_dt, pop_mfis, by = c("pop"))


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
        extracted_cd3 <- flowWorkspace::gh_pop_get_data(applied_gates, gatename)
        extracted_cd3 <- flowWorkspace::realize_view(extracted_cd3)
        extracted_cd3_ff <- flowWorkspace::cytoframe_to_flowFrame(extracted_cd3)

        # If an output directory is specified, write the extracted data to a new FCS file
        if (!is.na(outdir)) {
            filename_new <- sub(dir_fcs_compensated, "", fcs_x)
            # Rename the file properly
            filename_new <- sub("_ungated_", filename_ungated_into, filename_new)
            path_new <- file.path(outdir, filename_new)
            dir.create(dirname(path_new), recursive = TRUE, showWarnings = FALSE)
            flowCore::write.FCS(extracted_cd3_ff, path_new)

            if (verbose) {
                cat("   Wrote ", path_new, "\n")
            }
        }
    }

    # Remove populations that are not part of the usual gating strategy
    counts_complete <- counts_complete |> dplyr::filter(
        !pop %in% exclude_pops
    )
    # If an output directory is specified, write the counts to a CSV file
    if (!is.na(outdir)) {
        counts_complete |> data.table::fwrite(file.path(outdir, "counts_percent_MFI.csv"))
    }
    # Return the complete counts
    return(counts_complete)
}
