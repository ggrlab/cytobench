
#' Export fcs files
#'
#' Export a list of matrices, usually ``data.table``s.
#'
#' @param matrix_list
#'  List of matrices, should all have the same number of columns.
#' @param safety_scaling
#' An .fcs file is generated from the maximum cells (high and low). This safety scaling
#' is soleily for the range of the .fcs files. This does not change any cell values.
#' @param safety_shift
#' CAREFUL! Introducing a value here is changing all cell values. It adds this value
#' to all cells. This is "necessary" as Kaluza cannot visualize negative values on a
#' linear scale in .fcs files.
#' @param outdir
#' Where should the list of matrices be written to. Filenames come from names(matrix_list).
#' @param use.names Should the column names of each matrix be used when concattenating
#' the extreme .fcs?
#' @param verbose Verbose output?
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
#'
#' @param extreme_template
#' If you have an extreme sample which you want to use as a template for all other samples,
#' you can provide it here. Usually, this would be the returned fcs_extreme_copy from a
#' previous run of this function.
#' Essentially, it needs at least 2 cells (high and low) as "extreme" values.
#' The cells are replaced by the actual cells given in the matrix_list.
#'
#' @param cytname
#' The name of the cytometer. This is a custom keyword that we use to track the source. It is e.g. used by Kaluza to identify the default parameters for the cytometer.
#' @return fcs_extreme_copy:
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
                if (verbose) {
                    cat(".")
                }
                read_files_bound <- rbind(read_files_bound, x, use.names = FALSE)
            }
        }

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

        .SD <- NULL #  only for linting
        extreme_template <- rbind(
            read_files_bound[, lapply(.SD, max)],
            read_files_bound[, lapply(.SD, min)]
        )
        # have a security net of 20% plus 100 such that no values are negative (for Kaluza)
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
        flowCore::keyword(fcs_extreme)[["CODE_SOURCE"]] <- "2023-04-11_SSS_04/src/_functions/export_fcs.R"

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
