#' Read Custom Spillover Matrix and Autofluorescence Vector from FCS File
#'
#' Be careful when using it, this is heavily dependent on the FCS file metadata structure.
#'
#' This helper function extracts a custom spillover matrix and autofluorescence vector
#' from an FCS file or a `flowFrame` object. It is primarily intended for use with
#' `cmv_helper_compensation()` and related CMV data processing.
#'
#' The function expects the spillover matrix and autofluorescence vector to be encoded in
#' custom keywords within the FCS metadata, following a convention like
#' `"spillover.manual"` and `"spillover_autofluorescence.manual"`.
#'
#' @param fcs Either a character path to an FCS file or a `flowFrame` object.
#' @param custom_spillover_keyword Character. The keyword used to retrieve the custom
#'   spillover matrix and autofluorescence vector. Defaults to `"spillover.manual"`.
#' @param original_spillover_keyword Character. The keyword used to identify the original
#'   spillover matrix for determining channel order. Defaults to `"$SPILLOVER"`.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{`spillover`}{A numeric matrix representing the spillover matrix.}
#'   \item{`autofluorescence_proportional`}{A named numeric vector with autofluorescence contributions per channel.}
#' }
#' @export
#' @keywords cytometry
#' @examples
#' \dontrun{
#' read.FCS_custom_spillover(
#'     fcs = outfile,
#'     custom_spillover_keyword = comp_x,
#'     original_spillover_keyword = spillover_keyword
#' )
#' }
read.FCS_custom_spillover <- function(fcs,
                                      custom_spillover_keyword = "spillover.manual",
                                      original_spillover_keyword = "$SPILLOVER") {
    if (is.character(fcs)) {
        # If input is a file path, read only the FCS header to access metadata
        fs_keywords <- flowCore::read.FCSheader(fcs)[[1]]

        # Convert the original spillover matrix from string to matrix form
        fs_keywords[original_spillover_keyword] <- list(
            string_to_spill(fs_keywords[[original_spillover_keyword]])
        )
    } else {
        # If input is a flowFrame, extract metadata directly
        fs_keywords <- flowCore::keyword(fcs)
    }

    # Ensure the custom keyword is properly prefixed
    if (!startsWith(custom_spillover_keyword, "spillover.")) {
        custom_spillover_keyword <- paste0("spillover.", custom_spillover_keyword)
    }

    # flowCore uses the following function to read the spillover matrix:
    # read_spillover <- string_to_spill(fs_keywords[[original_spillover_keyword]])
    spillovermat <- string_to_spill(fs_keywords[[custom_spillover_keyword]])

    # Handle autofluorescence: zero vector for original, or read from paired keyword
    if (custom_spillover_keyword == "spillover.original") {
        autofluorescence_proportional <- spillovermat[1, ]
        autofluorescence_proportional[TRUE] <- 0
    } else {
        autofluorescence_proportional <- as.numeric(
            utils::read.table(text = fs_keywords[[sub(
                "spillover",
                "spillover_autofluorescence",
                custom_spillover_keyword
            )]])
        )
        names(autofluorescence_proportional) <- colnames(spillovermat)
    }

    # Ensure the output spillover matrix matches the ordering of the original
    current_spillover <- fs_keywords[[original_spillover_keyword]]
    # Order the spillover matrix according to the current spillover matrix
    spillover_reordered <- spillovermat
    rownames(spillover_reordered) <- colnames(spillover_reordered)
    spillover_reordered <- spillover_reordered[colnames(current_spillover), colnames(current_spillover)]

    return(
        list(
            spillover = spillover_reordered,
            autofluorescence_proportional = autofluorescence_proportional[colnames(current_spillover)]
        )
    )
}
