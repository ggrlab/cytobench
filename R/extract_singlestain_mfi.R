
extract_singlestain_mfi <- function(fcs_dir = "data-raw/s001",
                                    regex_singlestain = "-(CD3-.*)|(none)\\.fcs",
                                    transform = function(x) {
                                        asinh(x / 1e3)
                                    },
                                    gating_set_file = NULL,
                                    gate_extract = NULL) {
    dir_files <- list.files(fcs_dir, full.names = TRUE, pattern = regex_singlestain)
    loaded_fcs <- sapply(dir_files, flowWorkspace::load_cytoset_from_fcs, simplify = FALSE)
    if (!all(is.null(gating_set_file))) {
        if ("character" %in% class(gating_set_file)) {
            gatingset <- flowWorkspace::load_gs(gating_set_file)
        } else {
            gatingset <- gating_set_file
        }
        # Then gate the fcs files
        loaded_fcs <- lapply(loaded_fcs, function(ff_x) {
            suppressMessages(
                ff_x_gated <- flowWorkspace::gh_apply_to_cs(gatingset, ff_x, compensation_source = "none")
            )
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
            gated_ff <- flowWorkspace::gs_pop_get_data(ff_x_gated, gate_extract) |>
                flowWorkspace::realize_view()
            return(gated_ff)
        })
    }
    relevant_mfis <- sapply(names(loaded_fcs), function(f_x) {
        ff_x <- flowWorkspace::cytoframe_to_flowFrame(loaded_fcs[[f_x]][[1]])
        nonempty_channel <- which(flowCore::markernames(ff_x) != "empty")
        if (length(nonempty_channel) > 1) {
            stop("More than one non-empty channel in ", f_x)
        }
        if (length(nonempty_channel) == 0) {
            # Then it is the unstained sample --> return the MFI of all channels
            mfis <- apply(flowCore::exprs(ff_x), 2, median)
            mfis_tib <- tibble::tibble(
                "feature" = names(mfis),
                "unstained" = mfis,
            )
        } else {
            # first cluster the relevant channel into TWO populations,
            values_nonempty <- flowCore::exprs(ff_x)[, names(nonempty_channel)]
            # then return the median of both populations for that channel
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

    single_stainings <- do.call(rbind, relevant_mfis[!sapply(relevant_mfis, function(x) nrow(x) > 1)])
    relevant_unstained <- relevant_mfis[sapply(relevant_mfis, function(x) nrow(x) > 1)]
    if (length(relevant_unstained) == 1) {
        # If there was only one unstained sample, everything was fine, otherwise we just report everything
        # but that was probably a problem
        names(relevant_unstained) <- "unstained"
    } else {
        warning("More than one unstained sample found, returning each MFI by filename")
    }

    all_unstained_joint <- Reduce(dplyr::left_join, relevant_unstained)

    return(dplyr::left_join(
        single_stainings,
        all_unstained_joint,
        by = "feature"
    ))
}
