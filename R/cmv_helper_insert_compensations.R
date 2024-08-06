#' Insert multiple compensations into a fcs file
#'
#' This function is a helper function to INSERT the compensation into the fcs files of the CMV data from
#' Glehr et al. [TO BE SUBMITTED].
#' @param fcs_filename
#' The path to the fcs file to insert the compensations into
#' @param named_compensation_files
#' A named list of files containing the compensation files. The names are the compensation types.
#' @param outdir
#' The directory to save the fcs files with the inserted compensations. If NA, the function returns the flowFrame object with the inserted compensations.
#' @param sub_filename
#' A vector of two strings to substitute in the fcs_filename. The first string is substituted by the second string. The final files are saved in \code{file.path(outdir, sub(sub_filename[1], sub_filename[2], fcs_filename))}
#' @param verbose
#' Whether to print out the path of the saved fcs file
#' @export
cmv_helper_insert_compensations <- function(fcs_filename,
                                            named_compensation_files = list(
                                                # 1. Manual compensation, separately for single-stain and panel samples
                                                # The manual compensation used  "DCtubes_auto" compensation as basis
                                                "manual" = NA,
                                                # 2. Autocompensation, calculated on regular Duraclone Compensation tubes
                                                "DCtubes_auto" = NA,
                                                # 3. Autocompensation, calculated on a sample-by-sample basis using single-stain samples
                                                "singlestain_auto" = NA,
                                                # 4. Manual compensation, done on a sample-by-sample basis using concattenated single-stain samples
                                                # The manual compensation used "singlestain_auto" compensation as basis
                                                "singlestain_manual" = NA
                                            ), outdir = file.path("res", "CT_UnappliedCompensations"),
                                            sub_filename = c("intermediate/CT", ""),
                                            verbose = TRUE) {
    ff <- flowCore::read.FCS(fcs_filename)
    original_spillover <- flowCore::keyword(ff)[["$SPILLOVER"]]
    loaded_comps <- lapply(named_compensation_files, function(x) {
        if (is.null(x)) {
            return(list(
                "spillover" = NULL,
                "autofluorescence_proportion" = NULL
            ))
        }
        if (!file.exists(x)) {
            stop("Compensation directory ", x, " not found for file ", fcs_filename, "\n")
        }
        compfiles <- list.files(x, full.names = TRUE)
        file_spillover <- compfiles[grepl("spillover", compfiles)]
        file_autofluorescence <- compfiles[grepl("autofluorescence", compfiles)]

        spillover_mat <- as.matrix(data.table::fread(file_spillover))
        # rownames(spillover_mat) <- colnames(spillover_mat)
        colnames(spillover_mat) <- paste0(colnames(spillover_mat), "-A")
        spillover_mat <- spillover_mat[, colnames(original_spillover)]

        autofluorescence_vec <- as.numeric(data.table::fread(file_autofluorescence))
        names(autofluorescence_vec) <- colnames(spillover_mat)
        autofluorescence_vec[colnames(original_spillover)]

        return(
            list(
                "spillover" = spillover_mat,
                "autofluorescence_proportion" = autofluorescence_vec
            )
        )
    })

    # 1. Save the original spillover matrix - probably coming from the machine itself
    flowCore::keyword(ff)[["spillover.original"]] <- original_spillover

    # 2. Save all the calculated spillover matrices and autofluorescence vectors
    for (compensation_keyword in names(loaded_comps)) {
        flowCore::keyword(ff)[[paste0("spillover.", compensation_keyword)]] <- loaded_comps[[compensation_keyword]]$spillover
        flowCore::keyword(ff)[[paste0("spillover_autofluorescence.", compensation_keyword)]] <- loaded_comps[[compensation_keyword]]$autofluorescence_proportion
    }

    # NOTE: Kaluza CANNOT load the autofluorescence from the fcs file.
    # So the compensation will be slightly different to OUR Kaluza compensation
    # with this.
    flowCore::keyword(ff)[["$SPILLOVER"]] <- loaded_comps$manual$spillover

    if (!is.na(outdir)) {
        if (!all(is.na(sub_filename))) {
            fcs_filename <- sub(sub_filename[1], sub_filename[2], fcs_filename)
        }
        outfile <- file.path(outdir, fcs_filename)
        dir.create(
            dirname(outfile),
            recursive = TRUE,
            showWarnings = FALSE
        )
        flowCore::write.FCS(ff, outfile)
        if (verbose) {
            cat("FCS file containing ", paste0(names(named_compensation_files), collapse = ", "), " compensations saved to ", outfile, "\n")
        }
        # Test that reading the spillover+autofluorescence matrices works
        # This is how you would extract the compensations from the FCS file
        for (comp_x in c("spillover.original", names(loaded_comps))) {
            # cytobench::read.FCS_custom_spillover(outfile, comp_x)
            read.FCS_custom_spillover(outfile, comp_x)
        }
    } else {
        return(ff)
    }
}
