#' Gate Cells Using Predefined Gating Sets
#'
#' This function applies existing `GatingSet` objects to a `flowSet` (usually of compensated FCS files),
#' extracts gated populations (e.g., CD3+ cells), and optionally computes population statistics
#' such as counts, frequencies, and median fluorescence intensities (MFIs).
#'
#' @param flowset A `flowSet` or list of `flowFrame` objects representing compensated samples.
#'        These FCS files are expected to contain the compensation matrix, but not have compensation applied.
#' @param gatingset A named list of `GatingSet` objects, one per sample, or a single global `GatingSet`.
#' @param gatename Character or numeric. Name or index of the population to extract.
#'        Default is `"/Singlets/CD45+/CD3+"`. If numeric, the nth gate will be selected.
#' @param verbose Logical. Print progress information? Default: `TRUE`.
#' @param inplace Logical. If `TRUE`, modifies the `flowset` in place. Otherwise, a new `flowSet` is created.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{`counts`}{A `data.table` of per-population counts, frequencies, and MFIs across all samples.}
#'   \item{`flowset_gated`}{The updated `flowSet`, with each sample replaced by its gated `flowFrame`.}
#' }
#'
#' @details
#' - Assumes FCS filenames and gating set names match by sample ID (`sampleNames()`).
#' - If only one gating set is given, it is applied to all samples (useful for batch gating).
#' - Compensation is NOT applied again, assuming it's already been handled.
#' - All columns are temporarily treated as marker names to enable MFI extraction across the full feature set.
#'
#' @export
gate_cells <- function(flowset,
                       gatingset,
                       gatename = "/Singlets/CD45+/CD3+",
                       verbose = TRUE,
                       inplace = FALSE) {
    global_gating <- NULL
    if (length(gatingset) <= 1) {
        global_gating <- gatingset
    }

    # Ensure sample IDs in gating list match those in flowset
    if (is.null(global_gating)) {
        # Check if the panel files and gating files match
        ### Get the sample IDs from the directory names of the panel files
        fcs_panel_ids <- flowCore::sampleNames(flowset)
        gates_panel_ids <- names(gatingset)

        if (anyDuplicated(gates_panel_ids) != 0) {
            stop("There are duplicated sample IDs in the gating files, cannot assign them uniquely to the to-be-gated FCS files.")
        }
        ### Check if all panel files have a corresponding gating file
        if (!all(fcs_panel_ids %in% gates_panel_ids)) {
            # In the case of R1 and R2 data, each sample was measured on NAVIOS and Cytoflex,
            # thus it could be that the gating file of NAVIOS is missing.
            # However, this would then stop with an error in the following for-loop
            stop("Not all panel files have a corresponding gating file")
        }
    }

    # Work on a copy of the flowset unless inplace modification is requested
    if (!inplace) {
        flowset <- flowCore::flowSet(flowCore::flowSet_to_list(flowset))
    }

    counts_complete <- NA

    # Process each (compensated) FCS sample
    for (fcs_i in seq_along(flowset)) {
        fcs_x <- flowset[[fcs_i]]
        fcs_x_sID <- flowCore::sampleNames(flowset)[fcs_i]
        if (verbose) {
            cat("Processing", fcs_x_sID, "\n")
        }

        # Select appropriate gating set
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

        # Apply the gating set (compensation assumed applied already)
        # and suppress messages
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

        # Temporarily treat all columns as marker names for full MFI extraction
        # Set all columns as markernames such that pop.MFI applies to all of them
        all_cn_as_markernames <- c(
            stats::setNames(
                flowCore::colnames(applied_gates),
                flowCore::colnames(applied_gates)
            )
        )
        all_cn_as_markernames <- all_cn_as_markernames[all_cn_as_markernames != "empty"]
        all_cn_as_markernames <- all_cn_as_markernames[unique(names(all_cn_as_markernames))]
        # Extract the population MFIs
        tmp_gates <- flowWorkspace::gs_clone(applied_gates)
        flowCore::markernames(tmp_gates) <- all_cn_as_markernames

        # Get the counts and percentages for each population
        counts_freqs <- lapply(c("count", "percent"), function(statistic) {
            flowWorkspace::gs_pop_get_stats(
                tmp_gates,
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
        pop_mfis <- flowWorkspace::gh_pop_get_stats(tmp_gates, type = flowWorkspace::pop.MFI)
        counts_dt_MFI <- dplyr::left_join(counts_dt, pop_mfis, by = "pop")
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
        # Extract the gatename (e.g. CD3+) population data
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
