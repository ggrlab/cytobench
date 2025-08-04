#' Insert Multiple Compensations into an FCS File
#'
#' This helper function inserts multiple compensation matrices and autofluorescence vectors into
#' the metadata of an FCS file, following the structure used in the CMV dataset
#' (Glehr et al., [TO BE SUBMITTED]). These values can later be retrieved conveniently
#' with the `read.FCS_custom_spillover()` function. They can be directly _applied_ with
#' `cmv_helper_compensation()`.
#'
#' @param fcs_filename Character. Path to the FCS file into which compensations should be inserted.
#' @param named_compensation_files Named list of directories. Each directory must contain two files:
#'   a "spillover" matrix and an "autofluorescence" vector file. Files are assigned with regex.
#' The list names define the compensation types  (e.g., `"manual"`, `"DCtubes_auto"`).
#' Default provides placeholders for the CMV paper compensations. Empty values (`NULL`) indicate
#' that this compensation type should not be inserted.
#' @param outdir Character or NA. Directory to save the updated FCS file. If `NA`, the function
#'   returns the modified `flowFrame` instead of saving.
#' @param sub_filename Character vector of length 2. If provided, performs a `sub()` on the input filename
#'   to determine the output path: `file.path(outdir, sub(sub_filename[1], sub_filename[2], fcs_filename))`.
#'   May be NA.
#' @param spillover_keyword Character. Keyword used to locate the original spillover matrix within the FCS metadata.
#'   Defaults to `"$SPILLOVER"`. Necessary for ordering the new spillover matrices correctly and
#'   behaves as automatic check if all necessary columns of the new spillover matrices are present.
#'
#' The ordering of the spillover matrix must be exactly the same as the ordering of the markers in the expression matrix, otherwise KALUZA(!!) will apply the wrong compensation values.
#' @param verbose Logical. If `TRUE`, prints the output path and inserted compensation types.
#'
#' @return If `outdir` is `NA`, returns the updated `flowFrame`. Otherwise, writes the modified FCS file and returns invisible.
#' @export
cmv_helper_insert_compensations <- function(fcs_filename,
                                            named_compensation_files = list(
                                                # 1. Manual compensation, separately for single-stain and panel samples
                                                # The manual compensation used  "DCtubes_auto" compensation as basis
                                                "manual" = NULL,
                                                # 2. Autocompensation, calculated on regular Duraclone Compensation tubes
                                                "DCtubes_auto" = NULL,
                                                # 3. Autocompensation, calculated on a sample-by-sample basis using single-stain samples
                                                "singlestain_auto" = NULL,
                                                # 4. Manual compensation, done on a sample-by-sample basis using concattenated single-stain samples
                                                # The manual compensation used "singlestain_auto" compensation as basis
                                                "singlestain_manual" = NULL
                                            ), outdir = file.path("res", "CT_UnappliedCompensations"),
                                            sub_filename = c("intermediate/CT", ""),
                                            spillover_keyword = "$SPILLOVER",
                                            verbose = TRUE) {
    # Read original FCS file and extract original spillover matrix
    ff <- flowCore::read.FCS(fcs_filename)
    original_spillover <- flowCore::keyword(ff)[[spillover_keyword]]
    if (is.null(original_spillover)) {
        stop("Spillover matrix not found in FCS file: ", fcs_filename, "\n")
    }

    # Load compensation matrices and autofluorescence vectors from files
    loaded_comps <- lapply(named_compensation_files, function(x) {
        if (is.null(x)) {
            return(list(
                "spillover" = NULL,
                "autofluorescence_proportion" = NULL
            ))
        }
        if (!file.exists(x)) {
            stop("Compensation directory not found: ", x, " (for file ", fcs_filename, ")\n")
        }

        compfiles <- list.files(x, full.names = TRUE)
        file_spillover <- compfiles[grepl("spillover", compfiles)]
        file_autofluorescence <- compfiles[grepl("autofluorescence", compfiles)]

        # Load and format spillover matrix
        spillover_mat <- as.matrix(data.table::fread(file_spillover))
        colnames(spillover_mat) <- paste0(colnames(spillover_mat), "-A")
        rownames(spillover_mat) <- colnames(spillover_mat)
        spillover_mat <- spillover_mat[colnames(original_spillover), colnames(original_spillover)]

        # Load and format autofluorescence vector
        # Fread to have the column names
        autofluorescence_vec <- data.table::fread(file_autofluorescence)
        # Take first row with dropping (would be implicit)
        # to create a (named) vector)
        autofluorescence_vec <- as.matrix(autofluorescence_vec)[1, , drop = TRUE]
        names(autofluorescence_vec) <- paste0(names(autofluorescence_vec), "-A")
        # Then sort the autofluorescence vector to match the spillover matrix
        # which was sorted according to the original spillover matrix.
        autofluorescence_vec <- autofluorescence_vec[colnames(spillover_mat)]

        return(list(
            "spillover" = spillover_mat,
            "autofluorescence_proportion" = autofluorescence_vec
        ))
    })

    # Save the original spillover matrix as a separate metadata keyword
    flowCore::keyword(ff)[["spillover.original"]] <- flowCore:::spill2txt(original_spillover)

    # Insert all loaded compensations into metadata
    for (compensation_keyword in names(loaded_comps)) {
        flowCore::keyword(ff)[[paste0("spillover.", compensation_keyword)]] <-
            flowCore:::spill2txt(loaded_comps[[compensation_keyword]]$spillover)
        flowCore::keyword(ff)[[paste0("spillover_autofluorescence.", compensation_keyword)]] <-
            loaded_comps[[compensation_keyword]]$autofluorescence_proportion
    }

    # NOTE: Kaluza CANNOT load the autofluorescence from the fcs file.
    # So the compensation will be slightly different to OUR Kaluza compensation
    # with this.
    # The ordering of the spillover matrix must be exactly the same as the ordering of the markers in the expression matrix, otherwise KALUZA(!!) will apply the wrong compensation values.
    flowCore::keyword(ff)[[spillover_keyword]] <- loaded_comps$manual$spillover[
        names(flowCore::markernames(ff)),
        names(flowCore::markernames(ff))
    ]

    # Save to disk or return updated flowFrame
    if (!is.na(outdir)) {
        if (!all(is.na(sub_filename))) {
            fcs_filename <- sub(sub_filename[1], sub_filename[2], fcs_filename)
        }
        outfile <- file.path(outdir, fcs_filename)
        dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
        flowCore::write.FCS(ff, outfile)
        if (verbose) {
            cat(
                "FCS file containing ",
                paste0(names(named_compensation_files), collapse = ", "),
                " compensations saved to ",
                outfile, "\n"
            )
        }

        # Test that reading the spillover+autofluorescence matrices works
        # This is how you would extract the compensations from the FCS file
        for (comp_x in c("spillover.original", names(loaded_comps))) {
            # cytobench::read.FCS_custom_spillover(outfile, comp_x)
            read.FCS_custom_spillover(
                fcs = outfile,
                custom_spillover_keyword = comp_x,
                original_spillover_keyword = spillover_keyword
            )
        }
        return(invisible(ff))
    } else {
        return(ff)
    }
}
