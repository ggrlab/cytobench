
#' Export Flow Cytometry Data as FCS Files
#'
#' This function exports a list of matrices or `data.table`s as `.fcs` files using
#' `flowCore::flowFrame()` and `flowCore::write.FCS()`. It ensures FCS-compatible formatting,
#' adds optional safety scaling and shifting, handles marker names, and builds a shared
#' extreme template sample to define channel ranges for all files.
#' 
#' This is particularly useful for exporting flow cytometry data for use in software like Kaluza,
#' which reads the metadata to set up the plots. Especially the maximum and minimum channel
#' values are important for Kaluza to visualize the data correctly.
#'
#' @param matrix_list A named list of matrices or `data.table`s, all with the same number of columns.
#' @param safety_scaling Numeric. A scaling factor applied to the max/min range template. Does not modify actual data.
#'   Default is `1.20` (20% safety buffer).
#' @param safety_shift Numeric. CAREFUL - MODIFIES DATA AT SAVE! 
#'  An additive shift applied to all cell values before export. 
#'  Necessary to avoid negative values in linear scale visualizations in tools like Kaluza. Default is `0`.
#' @param outdir Character. Output directory for the `.fcs` files. Filenames are derived from `names(matrix_list)`.
#' @param use.names Logical. If `TRUE`, retains column names during the concatenation of extreme samples.
#' @param verbose Logical. Whether to print progress messages. Default is `TRUE`.
#' @param new_colnames_to In our usecase, there is often the case where we have:
#'
#' Sample1:
#'  FL1 CD3 FITC,FL2 CD3 PE,FL3 CD3 ECD,FL4 CD3 PC5.5,FL5 CD3 PC7,FL6 CD3 APC,FL7 CD3 AF700,FL8 CD3 AA750,FL9 CD3 PB,FL10 CD3 KrO
#' Sample2:
#'  FL1 CD3 FITC,FL2 empty,FL3 empty,FL4 empty,FL5 empty,FL6 empty,FL7 empty,FL8 empty,FL9 empty,FL10 empty
#'
#' In this case, we want to use the column names of Sample1, but we want to retain the
#' column names of Sample2. The original column names of each sample are saved in either
#' the description or the colnames of the .fcs file. ``Colnames`` are the relevant
#' within Kaluza where the marker names are extracted (and gates applied).
#'
#' @param new_colnames Can be single element numeric or character - referring to the
#' sample whose column names should be used.
#' Alternatively, give a named vector of strings with the same length as the matrix has columns, e.g.:
#'
#' c(
#'     "FSC-A" = NA,
#'     "SSC-A" = NA,
#'     "FL1-A" = "marker",
#'     "FL4-A" = "marker",
#'     "FL5-A" = "marker",
#'     "FL7-A" = "marker",
#'     "FL8-A" = "marker",
#'     "FL9-A" = "marker",
#'     "FL10-A" = "marker",
#'     "FL11-A" = "marker",
#'     "FL12-A" = "marker",
#'     "FL13-A" = "marker",
#'     "FL14-A" = "marker",
#'     "FL16-A" = "marker",
#'     "FL21-A" = "marker",
#'     "FSC-Width" = NA,
#'     "TIME" = NA
#' )
#'
#' Every value which is called "marker" will be used as a marker name setting
#' flowCore::markernames() IN ORDER! The other values are ignored.
#'
#' @param extreme_template A matrix or `flowFrame` to use as a shared max/min template across all samples.
#' 
#' If you have an extreme sample which you want to use as a template for all other samples,
#' you can provide it here. Usually, this would be the returned fcs_extreme_copy from a
#' previous run of this function.
#' Essentially, it needs at least 2 cells (high and low) as "extreme" values.
#' The cells are replaced by the actual cells given in the matrix_list.
#' @param cytname Character. 
#' The name of the cytometer. Stored in the `$CYT` keyword, used by Kaluza for default instrument settings.
#' You can use this to 1) track the source of the data and 2) have one cytname for each extreme template you used. 
#' 
#' @return A `flowFrame` (`fcs_extreme_copy`) representing the extreme template.
#' A copy of the extreme_template which was used to save all samples. The fcs-file
#' ranges for all samples come from this extreme sample.
#' 
#' @export
export_fcs <- function(matrix_list,
                       safety_scaling = 1.20,
                       safety_shift = 0,
                       outdir = ".",
                       use.names = FALSE,
                       verbose = TRUE,
                       new_colnames_to = c("description", "colnames"),
                       new_colnames = 1,
                       extreme_template = NULL,
                       cytname = "cytobench_exportedFCS") {

    matrix_list_dt <- lapply(matrix_list, data.table::as.data.table)
    
    if (verbose) {
        cat("Concattenating extreme values of all samples with data.table\n")
    }
    
    # Build extreme template if not provided
    if (is.null(extreme_template)) {
        matrix_list_dt_extremes <- lapply(matrix_list_dt, function(x) {
            rbind(
                x[, lapply(.SD, max)],
                x[, lapply(.SD, min)]
            )
        })

        if (use.names) {
            read_files_bound <- do.call(rbind, matrix_list_dt_extremes)
        } else {
            read_files_bound <- matrix_list_dt_extremes[[1]]
            for (x in matrix_list_dt_extremes[-1]) {
                if (verbose) cat(".")
                read_files_bound <- rbind(read_files_bound, x, use.names = FALSE)
            }
        }

        # Assign final column names for the extreme template
        if (length(new_colnames) == 1 &&
            (is.numeric(new_colnames[1]) || new_colnames %in% names(matrix_list_dt))) {
            # If new_colnames = 1, this is default behaviour with rbind(...., use.names=FALSE)
            colnames(read_files_bound) <- colnames(matrix_list_dt[[new_colnames]])
        } else {
            if (length(new_colnames) != ncol(read_files_bound)) {
                stop(
                    paste0(
                        "Must either give an integer or single string referring to the ",
                        "sample whose column names should be used. Alternatively, ",
                        "give a specific vector of strings with the same length as the",
                        "matrix has columns."
                    )
                )
            }
            colnames(read_files_bound) <- names(new_colnames)
        }

        # Create extreme matrix with safety scaling/shifting
        .SD <- NULL #  only for linting
        extreme_template <- rbind(
            read_files_bound[, lapply(.SD, max)],
            read_files_bound[, lapply(.SD, min)]
        )
        # E.g. have a security net of 20% plus 100 such that no values are negative (for Kaluza)
        extreme_template_safeties <- as.matrix(extreme_template) * safety_scaling + safety_shift
        fcs_extreme <- flowCore::flowFrame(extreme_template_safeties)
        fcs_extreme_copy <- flowCore::flowFrame(extreme_template_safeties)
    } else {
        if ("flowFrame" %in% class(extreme_template)) {
            extreme_template <- flowCore::exprs(extreme_template)
        }
        fcs_extreme <- flowCore::flowFrame(as.matrix(extreme_template) + safety_shift)
        fcs_extreme_copy <- flowCore::flowFrame(as.matrix(extreme_template) + safety_shift)
    }

    # Loop through each sample and export as .fcs
    for (file_x in names(matrix_list_dt)) {
        current_file <- file.path(outdir, paste0(file_x, ".fcs"))
        dir.create(dirname(current_file), recursive = TRUE, showWarnings = FALSE)

        # With this replacement of colnames, the "Description" will always be the names of
        # the extreme_template
        if (new_colnames_to[1] == "description") {
            if (length(new_colnames) == ncol(matrix_list_dt[[file_x]])) {
                flowCore::markernames(fcs_extreme)[TRUE] <- names(new_colnames)[new_colnames == "marker"]
            } else {
                flowCore::markernames(fcs_extreme)[TRUE] <- colnames(matrix_list_dt[[file_x]])[colnames(matrix_list_dt[[file_x]]) != "TIME"]
            }
            colnames(matrix_list_dt[[file_x]]) <- flowCore::colnames(fcs_extreme)
        } else {
            flowCore::colnames(fcs_extreme) <- colnames(matrix_list_dt[[file_x]])
        }

        # Set identifier names for the FCS file
        # $FIL is the official standard
        # GUID and FILENAME are not standard, but we (and others!) already use them
        # GUID seems to come from flowCore.
        flowCore::keyword(fcs_extreme)[["$FIL"]] <- file_x
        flowCore::keyword(fcs_extreme)[["GUID"]] <- file_x
        flowCore::keyword(fcs_extreme)[["FILENAME"]] <- file_x
        flowCore::keyword(fcs_extreme)[["$CYT"]] <- cytname
        # This is a custom (by us) keyword that we use to track
        # the source of how the fcs file was generated
        flowCore::keyword(fcs_extreme)[["CODE_SOURCE"]] <- paste0("cytobench_", utils::packageVersion("cytobench"))

        # Assign expression values and write FCS
        # Replace the data with the matrix
        Biobase::exprs(fcs_extreme) <- as.matrix(matrix_list_dt[[file_x]]) + safety_shift
        flowCore::write.FCS(fcs_extreme, current_file)

        if (verbose) {
            cat("\nWrote ", current_file)
        }
    }

    if (verbose) {
        cat("\n")
    }

    return(fcs_extreme_copy)
}
