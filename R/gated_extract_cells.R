gated_extract_cells <- function(gatingset, gatenames = c("root"), return_x = c("flowset", "cytoset", "gatinghierarchy")) {
    # Get the counts and percentages for each population
    counts_freqs <- lapply(c("count", "percent"), function(statistic) {
        flowWorkspace::gs_pop_get_stats(
            gatingset,
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

    # Combine the counts with the complete counts
    if (all(is.na(counts_complete))) {
        counts_complete <- counts_dt_MFI
    } else {
        counts_complete <- dplyr::bind_rows(counts_complete, counts_dt_MFI)
    }

    # If gatename is numeric, get the corresponding population path
    pop_paths <- flowWorkspace::gs_get_pop_paths(gatingset)
    gatenames <- sapply(gatenames, function(gatename) {
        if (is.numeric(gatename)) {
            pop_paths[[gatename]]
        } else {
            gatename
        }
    })


    # Extract the gatename (e.g. CD3+) population data
    if (return_x[1] == "gatinghierarchy") {
        return(gatingset)
    }
    final_extracted <- sapply(gatenames, simplify = FALSE, function(gatename) {
        extracted <- flowWorkspace::gh_pop_get_data(gatingset, gatename)
        extracted <- flowWorkspace::realize_view(extracted)
        if (return_x[1] == "cytoset") {
            extracted
        } else {
            flowWorkspace::cytoframe_to_flowFrame(extracted)
        }
    })
    return(final_extracted)
}
