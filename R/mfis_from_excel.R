
#' Extract MFIs from Excel file
#'
#'
#' @param excel_file
#' The excel file must have a very specific format with the following columns:
#'      EXP
#'      CD3 LOT
#'      MACHINE
#'      CHANNEL
#'      DYE
#'      Name of the data file containing the data set
#'      MFI FSC
#'      MFI SSC
#'      MFI POSITIVE POPULATION
#'      MFI NEGATIVE POPULATION (GG)
#'      MFI NEGATIVE POPULATION  UNSTAINED
#'
#' Every single staining sample should be delimited by exactly one blank row.
#' The two machines "MFI NAVIOS" and "MFI EX" should be named sheets.
#' @param positive_mfi_colname
#' Column name of the excel sheets with the positive MFI
#' @param negative_mfi_colname
#' Column name of the excel sheets with the negative MFI.
#' @return
#'  List ("intern", "extern") of matrices with MFIs per sample, looking like
#' this:
#' $intern
#' $intern$Align_01.Lot_A.NAVIOS
#'                              FL1 CD3 FITC FL2 CD3 PE FL3 CD3 ECD FL4 CD3 PC5.5 FL5 CD3 PC7 FL6 CD3 APC FL7 CD3 AF700 FL8 CD3 AA750 FL9 CD3 PB FL10 CD3 KrO
#' MFI POSITIVE POPULATION          37488.82   131657.6   139639.11     790976.38    451150.3    54247.44      37088.07      24695.75   27504.88     10440.54
#' MFI NEGATIVE POPULATION (GG)       219.82      298.0      437.83        615.24       364.9       66.68        164.55         98.35     824.71       792.06
#' @examples
#' # See tests
#' @export
mfis_from_excel <- function(excel_file = "Pre_Arcsinh_Median_FI.xlsx",
                            positive_mfi_colname = "MFI POSITIVE POPULATION",
                            negative_mfi_colname = c("extern" = "MFI NEGATIVE POPULATION  UNSTAINED")) {
    # read all sheets
    read_mfi_pre_asinh_persample <- read_excel_allsheets(
        # filename = file.path(DIR_raw_data_unified, "Pre_Arcsinh_Median_FI.xlsx")
        filename = excel_file
    )

    ## Split samples
    sample_mfis <- lapply(read_mfi_pre_asinh_persample, function(x) {
        x$is_blank <- apply(x, 1, function(y) all(is.na(y)))
        # FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE
        # 0  0  0  0  0  0  0  0  0  0  0  0  1  1
        x$sample_id <- cumsum(x$is_blank)
        # Remove (previously) blank rows
        x_noblank <- x[!x$is_blank, ]
        return(split(x_noblank, x_noblank$sample_id))
    })
    sample_mfis <- unlist(sample_mfis, recursive = FALSE)

    # dims <- vapply(sample_mfis, dim, FUN.VALUE = numeric(2))
    # if (!all(apply(dims, 2, function(x) x == c(12, 13)))) {
    #     stop(paste0(
    #         "Not all samples have 12 rows and 13 columns. Check the data.",
    #         "12 rows for 11 single-stained markers and one unstained control. "
    #     ))
    # }

    mfis_formatted <- list()
    if (length(negative_mfi_colname) > 1) {
        if (is.null(names(negative_mfi_colname))) {
            stop("If multiple negative MFI columns are given, they must be named.")
        }
    }
    for (negative_mfi_col_X in negative_mfi_colname) {
        tmp_formatted <- lapply(sample_mfis, function(x) {
            # Remove the last two rows:
            #   Last is the verify-sample
            #   Second to last is the unstained control
            x_stainings <- x[1:(nrow(x) - 2), ]
            mfi_mat <- t(as.matrix(
                x_stainings[, c(positive_mfi_colname, negative_mfi_col_X)]
            ))
            mfi_mat <- mfi_mat[, -c(11, 12)]
            colnames(mfi_mat) <- paste0(x_stainings$CHANNEL, " CD3 ", x_stainings$DYE)
            return(mfi_mat)
        })
        names(tmp_formatted) <- vapply(sample_mfis, function(x) {
            x <- x[1, , drop = FALSE]
            return(paste(x[["EXP"]], x[["CD3 LOT"]], x[["MACHINE"]], sep = "."))
        }, FUN.VALUE = character(1))
        mfis_formatted[[negative_mfi_col_X]] <- tmp_formatted
    }
    names(mfis_formatted) <- names(negative_mfi_colname)

    if (any(vapply(mfis_formatted, function(x) {
        length(x) != length(unique(names(x)))
    }, FUN.VALUE = logical(1)))) {
        stop("Not all samples have unique names. Check the data.")
    }

    if (length(mfis_formatted) == 1) {
        mfis_formatted <- mfis_formatted[[1]]
    }

    return(mfis_formatted)
}
