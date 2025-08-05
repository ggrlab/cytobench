#' Extract Median Fluorescence Intensities (MFIs) from Excel File
#'
#' This function reads an Excel file containing MFIs for single-stained cytometry samples
#' and returns a list of MFI matrices, grouped by sample and machine. Each matrix contains
#' the positive and negative MFIs for each channel-marker combination.
#'
#' @param excel_file Character. Path to the Excel file. The file must contain a very specific format
#' with at least the following columns:
#' \itemize{
#'   \item \code{EXP}
#'   \item \code{CD3 LOT}
#'   \item \code{MACHINE}
#'   \item \code{CHANNEL}
#'   \item \code{DYE}
#'   \item \code{Name of the data file containing the data set}
#'   \item \code{MFI FSC}
#'   \item \code{MFI SSC}
#'   \item \code{MFI POSITIVE POPULATION}
#'   \item \code{MFI NEGATIVE POPULATION (GG)}
#'   \item \code{MFI NEGATIVE POPULATION  UNSTAINED}
#' }
#' Each sample should be separated by exactly one blank row.
#' Named sheets such as "MFI NAVIOS" or "MFI EX" are assumed.
#'
#' @param positive_mfi_colname Character. Column name in the Excel file for the positive MFI values.
#' Default is \code{"MFI POSITIVE POPULATION"}.
#'
#' @param negative_mfi_colname Named character vector. Column name(s) in the Excel file for the negative MFI values.
#' The names are used as the output list names (e.g. "extern", "intern").
#' Default is \code{c("extern" = "MFI NEGATIVE POPULATION  UNSTAINED")}.
#'
#' @return A named list (e.g. "intern", "extern") of per-sample MFI matrices. Each element of the list is itself a
#' list where each entry is a matrix with:
#' \describe{
#'   \item{Rows}{Two rows: one for positive MFI and one for negative MFI.}
#'   \item{Columns}{Channel-marker names (e.g. \code{"FL1 CD3 FITC"}, \code{"FL2 CD3 PE"}, etc.).}
#' }
#'
#' If only one negative column is provided, the output is a single named list (not wrapped in another list).
#'
#'  List ("intern", "extern") of matrices with MFIs per sample, looking like
#' this:
#' $intern
#' $intern$Align_01.Lot_A.NAVIOS
#'                              FL1 CD3 FITC FL2 CD3 PE FL3 CD3 ECD FL4 CD3 PC5.5 FL5 CD3 PC7 FL6 CD3 APC FL7 CD3 AF700 FL8 CD3 AA750 FL9 CD3 PB FL10 CD3 KrO
#' MFI POSITIVE POPULATION          37488.82   131657.6   139639.11     790976.38    451150.3    54247.44      37088.07      24695.75   27504.88     10440.54
#' MFI NEGATIVE POPULATION (GG)       219.82      298.0      437.83        615.24       364.9       66.68        164.55         98.35     824.71       792.06
#'
#' @examples
#' # See tests
#' @export
mfis_from_excel <- function(excel_file = "Pre_Arcsinh_Median_FI.xlsx",
                            positive_mfi_colname = "MFI POSITIVE POPULATION",
                            negative_mfi_colname = c("extern" = "MFI NEGATIVE POPULATION  UNSTAINED")) {
    # Read all sheets of the Excel file
    read_mfi_pre_asinh_persample <- read_excel_allsheets(filename = excel_file)

    # Split each sheet into individual samples based on blank rows
    sample_mfis <- lapply(read_mfi_pre_asinh_persample, function(x) {
        x$is_blank <- apply(x, 1, function(y) all(is.na(y)))
        x$sample_id <- cumsum(x$is_blank)
        # Remove (previously) blank rows
        x_noblank <- x[!x$is_blank, ]
        split(x_noblank, x_noblank$sample_id)
    }) |>
        unlist(recursive = FALSE)

    # Prepare formatted output
    mfis_formatted <- list()

    if (length(negative_mfi_colname) > 1 && is.null(names(negative_mfi_colname))) {
        stop("If multiple negative MFI columns are given, they must be named.")
    }

    for (neg_col in negative_mfi_colname) {
        tmp_formatted <- lapply(sample_mfis, function(x) {
            # Remove last two rows: verify-sample and unstained control
            x_stainings <- x[1:(nrow(x) - 2), ]
            mfi_mat <- t(as.matrix(x_stainings[, c(positive_mfi_colname, neg_col)]))
            # Drop channel 11 and 12 (verify and unstained)
            mfi_mat <- mfi_mat[, -c(11, 12)]
            colnames(mfi_mat) <- paste0(x_stainings$CHANNEL, " CD3 ", x_stainings$DYE)
            return(mfi_mat)
        })

        # Name each matrix using EXP.LOT.MACHINE
        names(tmp_formatted) <- vapply(sample_mfis, function(x) {
            x <- x[1, , drop = FALSE]
            return(paste(x[["EXP"]], x[["CD3 LOT"]], x[["MACHINE"]], sep = "."))
        }, FUN.VALUE = character(1))

        mfis_formatted[[neg_col]] <- tmp_formatted
    }

    names(mfis_formatted) <- names(negative_mfi_colname)

    # Sanity check: ensure unique sample names
    if (any(vapply(mfis_formatted, function(x) {
        length(x) != length(unique(names(x)))
    }, FUN.VALUE = logical(1)))) {
        stop("Not all samples have unique names. Check the Excel file.")
    }

    # Unwrap if only one entry
    if (length(mfis_formatted) == 1) {
        mfis_formatted <- mfis_formatted[[1]]
    }

    return(mfis_formatted)
}
