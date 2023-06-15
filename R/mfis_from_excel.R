
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
#' Every single staining sample should be delimited by at least one blank row.
#' The two machines "MFI NAVIOS" and "MFI EX" should be named sheets.
#'
#' @return
#'  List ("intern", "extern") of matrices with MFIs per sample
#' @examples
mfis_from_excel <- function(excel_file = "Pre_Arcsinh_Median_FI.xlsx") {
    # read all sheets
    read_mfi_pre_asinh_persample <- read_excel_allsheets(
        # filename = file.path(DIR_raw_data_unified, "Pre_Arcsinh_Median_FI.xlsx")
        filename = excel_file
    )
    # exclude irrelevant columns
    read_mfi_pre_asinh_persample <- lapply(read_mfi_pre_asinh_persample, function(x) {
        x[, 1:which(colnames(x) == "MFI NEGATIVE POPULATION  UNSTAINED")]
    })
    # set proper names
    names(read_mfi_pre_asinh_persample) <- lapply(read_mfi_pre_asinh_persample, function(x) {
        return(paste0(x[["EXP"]][1], ifelse(x[["MACHINE"]][1] == "NAVIOS EX", "_EX", "")))
    })
    names(read_mfi_pre_asinh_persample) <- c(
        read_mfi_pre_asinh_persample[[1]][["EXP"]][1],
        paste0(read_mfi_pre_asinh_persample[[2]][["EXP"]][1], "_EX")
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

    dims <- vapply(sample_mfis, dim, FUN.VALUE = numeric(2))

    if (!all(apply(dims, 2, function(x) x == c(12, 13)))) {
        stop(paste0(
            "Not all samples have 12 rows and 13 columns. Check the data.",
            "12 rows for 11 single-stained markers and one unstained control. "
        ))
    }

    mfis_formatted <- list()
    for (negative_in_or_extern in c("intern", "extern")) {
        if (negative_in_or_extern == "intern") {
            negative_stainging_column <- "MFI NEGATIVE POPULATION (GG)"
        } else {
            negative_stainging_column <- "MFI NEGATIVE POPULATION  UNSTAINED"
        }
        tmp_formatted <- lapply(sample_mfis, function(x) {
            # Remove the last two rows:
            #   Last is the verify-sample
            #   Second to last is the unstained control
            x_stainings <- x[1:(nrow(x) - 2), ]
            mfi_mat <- t(as.matrix(
                x_stainings[, c("MFI POSITIVE POPULATION", negative_stainging_column)]
            ))
            mfi_mat <- mfi_mat[, -c(11, 12)]
            colnames(mfi_mat) <- paste0(x_stainings$CHANNEL, " CD3 ", x_stainings$DYE)
            return(mfi_mat)
        })
        names(tmp_formatted) <- vapply(sample_mfis, function(x) {
            x <- x[1, , drop = FALSE]
            return(paste(x[["EXP"]], x[["CD3 LOT"]], x[["MACHINE"]], sep = "."))
        }, FUN.VALUE = character(1))
        mfis_formatted[[negative_in_or_extern]] <- tmp_formatted
    }
    if (any(vapply(mfis_formatted, function(x) {
        length(x) != length(unique(names(x)))
    }, FUN.VALUE = logical(1)))) {
        stop("Not all samples have unique names. Check the data.")
    }
}
